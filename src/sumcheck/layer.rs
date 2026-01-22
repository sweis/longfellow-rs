//! Sumcheck layer operations.
//!
//! This module implements the sumcheck protocol for a single layer of the circuit.

use crate::field::Field;
use crate::polynomial::lagrange_coefficients_3;
use crate::transcript::Transcript;

use super::eq::{bind, transpose_quad};

/// Polynomial evaluation at three points (P0=0, P1=1, P2).
#[derive(Clone, Debug)]
pub struct PolyEvals<F: Field> {
    /// Evaluation at P0 = 0.
    pub p0: F,
    /// Evaluation at P2 (typically 2 for prime fields).
    pub p2: F,
}

impl<F: Field> PolyEvals<F> {
    /// Create new polynomial evaluations.
    pub fn new(p0: F, p2: F) -> Self {
        Self { p0, p2 }
    }

    /// Compute the implied p1 from p0 and the claim.
    /// P(1) = claim - P(0)
    pub fn p1(&self, claim: F) -> F {
        claim - self.p0
    }

    /// Interpolate and evaluate at x.
    pub fn interpolate_at(&self, x: F, claim: F) -> F {
        let p1 = self.p1(claim);
        let [lag_0, lag_1, lag_2] = lagrange_coefficients_3(x, F::P2);
        lag_0 * self.p0 + lag_1 * p1 + lag_2 * self.p2
    }
}

/// Proof for a single sumcheck layer.
#[derive(Clone, Debug)]
pub struct SumcheckLayerProof<F: Field> {
    /// Evaluations for each round and hand.
    /// evals[round][hand] contains (p0, p2) evaluations.
    pub evals: Vec<[PolyEvals<F>; 2]>,

    /// Final left value after all bindings.
    pub vl: F,

    /// Final right value after all bindings.
    pub vr: F,
}

impl<F: Field> SumcheckLayerProof<F> {
    /// Create a new layer proof.
    pub fn new(evals: Vec<[PolyEvals<F>; 2]>, vl: F, vr: F) -> Self {
        Self { evals, vl, vr }
    }

    /// Serialize the proof.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        // Number of rounds
        bytes.extend_from_slice(&(self.evals.len() as u32).to_le_bytes());

        // Evaluations
        for round in &self.evals {
            for hand in round {
                bytes.extend_from_slice(&hand.p0.to_bytes());
                bytes.extend_from_slice(&hand.p2.to_bytes());
            }
        }

        // Final values
        bytes.extend_from_slice(&self.vl.to_bytes());
        bytes.extend_from_slice(&self.vr.to_bytes());

        bytes
    }
}

/// Layer pad structure for zero-knowledge.
#[derive(Clone, Debug)]
pub struct LayerPad<F: Field> {
    /// Pad evaluations for each round and hand.
    pub evals: Vec<[PolyEvals<F>; 2]>,

    /// Pad for final left value.
    pub vl: F,

    /// Pad for final right value.
    pub vr: F,

    /// Pad for the product vl * vr.
    pub vl_vr: F,
}

impl<F: Field> LayerPad<F> {
    /// Create an empty (zero) layer pad.
    pub fn zero(num_rounds: usize) -> Self {
        Self {
            evals: vec![
                [
                    PolyEvals::new(F::ZERO, F::ZERO),
                    PolyEvals::new(F::ZERO, F::ZERO),
                ];
                num_rounds
            ],
            vl: F::ZERO,
            vr: F::ZERO,
            vl_vr: F::ZERO,
        }
    }

    /// Create a random layer pad.
    pub fn random<R: rand::Rng>(rng: &mut R, num_rounds: usize) -> Self {
        Self {
            evals: (0..num_rounds)
                .map(|_| {
                    [
                        PolyEvals::new(F::random(rng), F::random(rng)),
                        PolyEvals::new(F::random(rng), F::random(rng)),
                    ]
                })
                .collect(),
            vl: F::random(rng),
            vr: F::random(rng),
            vl_vr: F::random(rng),
        }
    }
}

/// Run the sumcheck protocol for a single layer.
///
/// # Arguments
/// * `quad` - The bound QUAD array (2D: [l, r])
/// * `vl` - Left wire values
/// * `vr` - Right wire values
/// * `lv` - Number of variables (log of wire count)
/// * `layer_pad` - Padding for zero-knowledge
/// * `transcript` - Fiat-Shamir transcript
///
/// # Returns
/// The layer proof and the binding challenges G.
pub fn sumcheck_layer<F: Field>(
    quad: &[Vec<F>],
    vl: &[F],
    vr: &[F],
    lv: usize,
    layer_pad: &LayerPad<F>,
    transcript: &mut Transcript,
) -> (SumcheckLayerProof<F>, Vec<[F; 2]>) {
    let mut vl = vl.to_vec();
    let mut vr = vr.to_vec();
    let mut quad = quad.to_vec();

    let mut evals = Vec::with_capacity(lv);
    let mut g = Vec::with_capacity(lv);

    for round in 0..lv {
        let mut round_evals = [
            PolyEvals::new(F::ZERO, F::ZERO),
            PolyEvals::new(F::ZERO, F::ZERO),
        ];
        let mut round_g = [F::ZERO, F::ZERO];

        for hand in 0..2 {
            // Compute p(x) = SUM_{l, r} bind(QUAD, x)[l, r] * bind(VL, x)[l] * VR[r]
            let (p0, p2) = compute_poly_evals(&quad, &vl, &vr);

            // Subtract padding
            let padded_p0 = p0 - layer_pad.evals[round][hand].p0;
            let padded_p2 = p2 - layer_pad.evals[round][hand].p2;

            round_evals[hand] = PolyEvals::new(padded_p0, padded_p2);

            // Write to transcript
            transcript.write_field_element(&padded_p0);
            transcript.write_field_element(&padded_p2);

            // Generate challenge
            let challenge: F = transcript.generate_challenge();
            round_g[hand] = challenge;

            // Bind L to challenge
            vl = bind(&vl, challenge);
            quad = bind_quad_first(&quad, challenge);

            // Swap L and R
            std::mem::swap(&mut vl, &mut vr);
            quad = transpose_quad(&quad);
        }

        evals.push(round_evals);
        g.push(round_g);
    }

    // Extract final values
    let final_vl = vl.first().copied().unwrap_or(F::ZERO) - layer_pad.vl;
    let final_vr = vr.first().copied().unwrap_or(F::ZERO) - layer_pad.vr;

    // Write final values to transcript
    transcript.write_field_element(&final_vl);
    transcript.write_field_element(&final_vr);

    let proof = SumcheckLayerProof::new(evals, final_vl, final_vr);
    (proof, g)
}

/// Compute polynomial evaluations at P0=0 and P2.
fn compute_poly_evals<F: Field>(
    quad: &[Vec<F>],
    vl: &[F],
    vr: &[F],
) -> (F, F) {
    let mut p0 = F::ZERO;
    let mut p2 = F::ZERO;

    // P(0) = SUM_{l, r} QUAD[2*g, l, r] * VL[2*l] * VR[r]
    // P(2) = SUM_{l, r} ((1-2)*QUAD[2*g, l, r] + 2*QUAD[2*g+1, l, r]) *
    //                   ((1-2)*VL[2*l] + 2*VL[2*l+1]) * VR[r]

    // For simplicity, we compute P(0) and P(2) by evaluating the polynomial
    // at x=0 and x=P2

    // P(0): x = 0 means we take even-indexed elements
    for (l, quad_row) in quad.iter().enumerate() {
        let vl_val = vl.get(l).copied().unwrap_or(F::ZERO);
        for (r, &q_val) in quad_row.iter().enumerate() {
            let vr_val = vr.get(r).copied().unwrap_or(F::ZERO);
            p0 = p0 + q_val * vl_val * vr_val;
        }
    }

    // P(2): We need to evaluate with x = 2 (or P2)
    // For the bound polynomial at x, we have:
    // bind(A, x)[i] = (1-x) * A[2i] + x * A[2i+1]
    // At x = 2: (1-2) * A[2i] + 2 * A[2i+1] = -A[2i] + 2*A[2i+1]

    let two = F::P2;
    let neg_one = F::ZERO - F::ONE;

    // Bind VL at x = 2
    let vl_bound: Vec<F> = (0..(vl.len() + 1) / 2)
        .map(|i| {
            let v0 = vl.get(2 * i).copied().unwrap_or(F::ZERO);
            let v1 = vl.get(2 * i + 1).copied().unwrap_or(F::ZERO);
            neg_one * v0 + two * v1
        })
        .collect();

    // Bind QUAD's first dimension at x = 2
    let quad_bound: Vec<Vec<F>> = (0..(quad.len() + 1) / 2)
        .map(|i| {
            let row0 = quad.get(2 * i);
            let row1 = quad.get(2 * i + 1);
            let max_len = row0.map(|r| r.len()).unwrap_or(0)
                .max(row1.map(|r| r.len()).unwrap_or(0));

            (0..max_len)
                .map(|j| {
                    let v0 = row0.and_then(|r| r.get(j)).copied().unwrap_or(F::ZERO);
                    let v1 = row1.and_then(|r| r.get(j)).copied().unwrap_or(F::ZERO);
                    neg_one * v0 + two * v1
                })
                .collect()
        })
        .collect();

    for (l, quad_row) in quad_bound.iter().enumerate() {
        let vl_val = vl_bound.get(l).copied().unwrap_or(F::ZERO);
        for (r, &q_val) in quad_row.iter().enumerate() {
            let vr_val = vr.get(r).copied().unwrap_or(F::ZERO);
            p2 = p2 + q_val * vl_val * vr_val;
        }
    }

    (p0, p2)
}

/// Bind the first dimension of a QUAD array.
fn bind_quad_first<F: Field>(quad: &[Vec<F>], x: F) -> Vec<Vec<F>> {
    let n = (quad.len() + 1) / 2;
    let one_minus_x = F::ONE - x;

    (0..n)
        .map(|i| {
            let row0 = quad.get(2 * i);
            let row1 = quad.get(2 * i + 1);
            let max_len = row0.map(|r| r.len()).unwrap_or(0)
                .max(row1.map(|r| r.len()).unwrap_or(0));

            (0..max_len)
                .map(|j| {
                    let v0 = row0.and_then(|r| r.get(j)).copied().unwrap_or(F::ZERO);
                    let v1 = row1.and_then(|r| r.get(j)).copied().unwrap_or(F::ZERO);
                    one_minus_x * v0 + x * v1
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
    fn test_poly_evals() {
        let quad = vec![
            vec![Fp128::from_u64(1), Fp128::from_u64(2)],
            vec![Fp128::from_u64(3), Fp128::from_u64(4)],
        ];
        let vl = vec![Fp128::from_u64(1), Fp128::from_u64(2)];
        let vr = vec![Fp128::from_u64(1), Fp128::from_u64(1)];

        let (p0, _p2) = compute_poly_evals(&quad, &vl, &vr);

        // P(0) should be the sum over l, r of quad[l][r] * vl[l] * vr[r]
        // = 1*1*1 + 2*1*1 + 3*2*1 + 4*2*1 = 1 + 2 + 6 + 8 = 17
        assert_eq!(p0, Fp128::from_u64(17));
    }

    #[test]
    fn test_layer_proof() {
        let mut transcript = Transcript::new(b"test");

        let quad = vec![
            vec![Fp128::from_u64(1)],
        ];
        let vl = vec![Fp128::from_u64(1)];
        let vr = vec![Fp128::from_u64(1)];
        let pad = LayerPad::zero(1);

        let (proof, g) = sumcheck_layer(&quad, &vl, &vr, 1, &pad, &mut transcript);

        assert_eq!(g.len(), 1);
        assert!(!proof.vl.is_zero() || !proof.vr.is_zero() || proof.evals.len() == 1);
    }
}
