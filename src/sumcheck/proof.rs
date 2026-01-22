//! Sumcheck circuit proof.
//!
//! This module implements the full sumcheck protocol over a layered circuit.

use crate::field::Field;
use crate::transcript::Transcript;

use super::layer::{sumcheck_layer, LayerPad, SumcheckLayerProof};

/// A layered circuit layer.
#[derive(Clone, Debug)]
pub struct CircuitLayer<F: Field> {
    /// The QUAD array for this layer: Q[g, l, r].
    /// Outer vec is g dimension, inner vecs are l, r dimensions.
    pub quad: Vec<Vec<Vec<F>>>,

    /// The Z array for in-circuit assertions (optional).
    pub z: Vec<Vec<Vec<F>>>,

    /// Log of the number of output wires.
    pub lv: usize,
}

impl<F: Field> CircuitLayer<F> {
    /// Create a new circuit layer.
    pub fn new(quad: Vec<Vec<Vec<F>>>, lv: usize) -> Self {
        Self {
            quad,
            z: vec![],
            lv,
        }
    }

    /// Create a layer with in-circuit assertions.
    pub fn with_assertions(quad: Vec<Vec<Vec<F>>>, z: Vec<Vec<Vec<F>>>, lv: usize) -> Self {
        Self { quad, z, lv }
    }

    /// Get the combined QZ = Q + beta * Z array.
    pub fn combined_quad(&self, beta: F) -> Vec<Vec<Vec<F>>> {
        if self.z.is_empty() {
            return self.quad.clone();
        }

        let mut qz = self.quad.clone();

        for (g, z_g) in self.z.iter().enumerate() {
            if g >= qz.len() {
                qz.resize_with(g + 1, Vec::new);
            }

            for (l, z_l) in z_g.iter().enumerate() {
                if l >= qz[g].len() {
                    qz[g].resize_with(l + 1, Vec::new);
                }

                for (r, &z_val) in z_l.iter().enumerate() {
                    if r >= qz[g][l].len() {
                        qz[g][l].resize(r + 1, F::ZERO);
                    }

                    qz[g][l][r] = qz[g][l][r] + beta * z_val;
                }
            }
        }

        qz
    }
}

/// Wire values for a layer.
#[derive(Clone, Debug)]
pub struct LayerWires<F: Field> {
    /// Left input wires.
    pub vl: Vec<F>,
    /// Right input wires.
    pub vr: Vec<F>,
}

impl<F: Field> LayerWires<F> {
    /// Create new layer wires.
    pub fn new(vl: Vec<F>, vr: Vec<F>) -> Self {
        Self { vl, vr }
    }

    /// Create from a single wire vector (vl = vr = v).
    pub fn from_single(v: Vec<F>) -> Self {
        Self {
            vl: v.clone(),
            vr: v,
        }
    }
}

/// Full sumcheck proof over a circuit.
#[derive(Clone, Debug)]
pub struct SumcheckProof<F: Field> {
    /// Layer proofs.
    pub layer_proofs: Vec<SumcheckLayerProof<F>>,
}

impl<F: Field> SumcheckProof<F> {
    /// Create a new sumcheck proof.
    pub fn new(layer_proofs: Vec<SumcheckLayerProof<F>>) -> Self {
        Self { layer_proofs }
    }

    /// Serialize the proof.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        bytes.extend_from_slice(&(self.layer_proofs.len() as u32).to_le_bytes());

        for proof in &self.layer_proofs {
            let layer_bytes = proof.to_bytes();
            bytes.extend_from_slice(&(layer_bytes.len() as u32).to_le_bytes());
            bytes.extend_from_slice(&layer_bytes);
        }

        bytes
    }
}

/// Padding for the entire circuit proof.
#[derive(Clone, Debug)]
pub struct CircuitPad<F: Field> {
    /// Layer pads.
    pub layer_pads: Vec<LayerPad<F>>,
}

impl<F: Field> CircuitPad<F> {
    /// Create zero padding.
    pub fn zero(layers: &[CircuitLayer<F>]) -> Self {
        Self {
            layer_pads: layers.iter().map(|l| LayerPad::zero(l.lv)).collect(),
        }
    }

    /// Create random padding.
    pub fn random<R: rand::Rng>(rng: &mut R, layers: &[CircuitLayer<F>]) -> Self {
        Self {
            layer_pads: layers.iter().map(|l| LayerPad::random(rng, l.lv)).collect(),
        }
    }
}

/// Run the sumcheck protocol over an entire circuit.
///
/// # Arguments
/// * `layers` - The circuit layers
/// * `wires` - Wire values for each layer
/// * `pad` - Circuit padding for zero-knowledge
/// * `transcript` - Fiat-Shamir transcript
///
/// # Returns
/// The sumcheck proof.
pub fn sumcheck_circuit<F: Field>(
    layers: &[CircuitLayer<F>],
    wires: &[LayerWires<F>],
    pad: &CircuitPad<F>,
    transcript: &mut Transcript,
) -> SumcheckProof<F> {
    let nl = layers.len();
    assert_eq!(wires.len(), nl);
    assert_eq!(pad.layer_pads.len(), nl);

    // Generate initial binding challenges
    let lv0 = if nl > 0 { layers[0].lv } else { 0 };
    let mut g: Vec<F> = transcript.generate_challenges(lv0);

    let mut layer_proofs = Vec::with_capacity(nl);

    for j in 0..nl {
        let layer = &layers[j];
        let layer_wires = &wires[j];
        let layer_pad = &pad.layer_pads[j];

        // Generate challenges
        let alpha: F = transcript.generate_challenge();
        let beta: F = transcript.generate_challenge();

        // Combine Q and Z
        let qz = layer.combined_quad(beta);

        // Bind the g dimension
        let bound_qz_0 = bind_quad_3d(&qz, &g);
        let bound_qz_1 = bind_quad_3d(&qz, &g);

        // Combine the two bindings
        let quad_2d = combine_bindings(&bound_qz_0, &bound_qz_1, alpha);

        // Run layer sumcheck
        let (proof, new_g) = sumcheck_layer(
            &quad_2d,
            &layer_wires.vl,
            &layer_wires.vr,
            layer.lv,
            layer_pad,
            transcript,
        );

        layer_proofs.push(proof);

        // Update g for next layer
        g = new_g.iter().flat_map(|arr| arr.iter().copied()).collect();
    }

    SumcheckProof::new(layer_proofs)
}

/// Bind a 3D QUAD array.
fn bind_quad_3d<F: Field>(quad: &[Vec<Vec<F>>], x: &[F]) -> Vec<Vec<F>> {
    let mut result: Vec<Vec<Vec<F>>> = quad.to_vec();

    for &xi in x {
        result = bind_3d_first(&result, xi);
    }

    // After binding all g variables, we should have a 2D array
    if result.len() == 1 {
        result.into_iter().next().unwrap_or_default()
    } else {
        // If not fully bound, flatten
        result.into_iter().flatten().collect()
    }
}

/// Bind the first dimension of a 3D array.
fn bind_3d_first<F: Field>(q: &[Vec<Vec<F>>], x: F) -> Vec<Vec<Vec<F>>> {
    let n = (q.len() + 1) / 2;
    let one_minus_x = F::ONE - x;

    (0..n)
        .map(|i| {
            let q0 = q.get(2 * i);
            let q1 = q.get(2 * i + 1);

            let l_size = q0.map(|v| v.len()).unwrap_or(0).max(q1.map(|v| v.len()).unwrap_or(0));

            (0..l_size)
                .map(|l| {
                    let row0 = q0.and_then(|v| v.get(l));
                    let row1 = q1.and_then(|v| v.get(l));

                    let r_size = row0.map(|r| r.len()).unwrap_or(0)
                        .max(row1.map(|r| r.len()).unwrap_or(0));

                    (0..r_size)
                        .map(|r| {
                            let v0 = row0.and_then(|v| v.get(r)).copied().unwrap_or(F::ZERO);
                            let v1 = row1.and_then(|v| v.get(r)).copied().unwrap_or(F::ZERO);
                            one_minus_x * v0 + x * v1
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}

/// Combine two bound QUAD arrays: result = a + alpha * b.
fn combine_bindings<F: Field>(
    a: &[Vec<F>],
    b: &[Vec<F>],
    alpha: F,
) -> Vec<Vec<F>> {
    let l_size = a.len().max(b.len());

    (0..l_size)
        .map(|l| {
            let row_a = a.get(l);
            let row_b = b.get(l);

            let r_size = row_a.map(|r| r.len()).unwrap_or(0)
                .max(row_b.map(|r| r.len()).unwrap_or(0));

            (0..r_size)
                .map(|r| {
                    let va = row_a.and_then(|row| row.get(r)).copied().unwrap_or(F::ZERO);
                    let vb = row_b.and_then(|row| row.get(r)).copied().unwrap_or(F::ZERO);
                    va + alpha * vb
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp128;

    #[test]
    fn test_circuit_layer() {
        // Simple identity layer: V[0] = V'[0] * V'[0]
        let quad = vec![vec![vec![Fp128::ONE]]];
        let layer = CircuitLayer::new(quad, 1);

        assert_eq!(layer.lv, 1);
    }

    #[test]
    fn test_sumcheck_simple() {
        let mut transcript = Transcript::new(b"test");

        // Simple one-layer circuit
        let quad = vec![vec![vec![Fp128::ONE]]];
        let layer = CircuitLayer::new(quad, 1);
        let layers = vec![layer];

        let wires = vec![LayerWires::from_single(vec![Fp128::from_u64(2)])];
        let pad = CircuitPad::zero(&layers);

        let proof = sumcheck_circuit(&layers, &wires, &pad, &mut transcript);

        assert_eq!(proof.layer_proofs.len(), 1);
    }
}
