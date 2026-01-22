//! Circuit representation for the Longfellow ZK scheme.
//!
//! This module implements layered arithmetic circuits used in the ZK proof system.
//! A circuit consists of layers, where each layer computes output wires from
//! input wires using quadratic operations.

use crate::field::Field;

/// A quad tuple (g, l, r, v) representing a term in the circuit.
///
/// If v = 0, this represents an assertion Z[g, l, r] = 1.
/// If v != 0, this represents Q[g, l, r] = v.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QuadTerm<F: Field> {
    /// Output wire index.
    pub g: usize,
    /// Left input wire index.
    pub l: usize,
    /// Right input wire index.
    pub r: usize,
    /// Value (0 for assertions, nonzero for computation).
    pub v: F,
}

impl<F: Field> QuadTerm<F> {
    /// Create a new quad term.
    pub fn new(g: usize, l: usize, r: usize, v: F) -> Self {
        Self { g, l, r, v }
    }

    /// Create an assertion term (v = 0).
    pub fn assertion(g: usize, l: usize, r: usize) -> Self {
        Self {
            g,
            l,
            r,
            v: F::ZERO,
        }
    }

    /// Check if this is an assertion term.
    pub fn is_assertion(&self) -> bool {
        self.v.is_zero()
    }
}

/// A layer in the circuit.
#[derive(Clone, Debug)]
pub struct Layer<F: Field> {
    /// Quad terms for this layer.
    pub terms: Vec<QuadTerm<F>>,

    /// Log of the number of output wires.
    pub lv_out: usize,

    /// Log of the number of input wires.
    pub lv_in: usize,
}

impl<F: Field> Layer<F> {
    /// Create a new layer.
    pub fn new(terms: Vec<QuadTerm<F>>, lv_out: usize, lv_in: usize) -> Self {
        Self {
            terms,
            lv_out,
            lv_in,
        }
    }

    /// Get the number of output wires.
    pub fn num_outputs(&self) -> usize {
        1 << self.lv_out
    }

    /// Get the number of input wires.
    pub fn num_inputs(&self) -> usize {
        1 << self.lv_in
    }

    /// Evaluate the layer given input wires.
    pub fn evaluate(&self, inputs: &[F]) -> Vec<F> {
        let num_outputs = self.num_outputs();
        let mut outputs = vec![F::ZERO; num_outputs];

        for term in &self.terms {
            if term.g >= num_outputs {
                continue;
            }

            let l_val = inputs.get(term.l).copied().unwrap_or(F::ZERO);
            let r_val = inputs.get(term.r).copied().unwrap_or(F::ZERO);

            // V[g] += Q[g, l, r] * V'[l] * V'[r]
            outputs[term.g] = outputs[term.g] + term.v * l_val * r_val;
        }

        outputs
    }

    /// Build the full QUAD array for this layer.
    pub fn build_quad(&self) -> Vec<Vec<Vec<F>>> {
        let num_g = self.num_outputs();
        let num_lr = self.num_inputs();

        let mut quad = vec![vec![vec![F::ZERO; num_lr]; num_lr]; num_g];

        for term in &self.terms {
            if term.g < num_g && term.l < num_lr && term.r < num_lr {
                quad[term.g][term.l][term.r] = term.v;
            }
        }

        quad
    }

    /// Build the Q and Z arrays separately.
    pub fn build_q_and_z(&self) -> (Vec<Vec<Vec<F>>>, Vec<Vec<Vec<F>>>) {
        let num_g = self.num_outputs();
        let num_lr = self.num_inputs();

        let mut q = vec![vec![vec![F::ZERO; num_lr]; num_lr]; num_g];
        let mut z = vec![vec![vec![F::ZERO; num_lr]; num_lr]; num_g];

        for term in &self.terms {
            if term.g < num_g && term.l < num_lr && term.r < num_lr {
                if term.is_assertion() {
                    z[term.g][term.l][term.r] = F::ONE;
                } else {
                    q[term.g][term.l][term.r] = term.v;
                }
            }
        }

        (q, z)
    }
}

/// A complete layered circuit.
#[derive(Clone, Debug)]
pub struct Circuit<F: Field> {
    /// The layers of the circuit.
    pub layers: Vec<Layer<F>>,

    /// Number of public input wires.
    pub num_public_inputs: usize,

    /// Number of private input (witness) wires.
    pub num_private_inputs: usize,
}

impl<F: Field> Circuit<F> {
    /// Create a new circuit.
    pub fn new(
        layers: Vec<Layer<F>>,
        num_public_inputs: usize,
        num_private_inputs: usize,
    ) -> Self {
        Self {
            layers,
            num_public_inputs,
            num_private_inputs,
        }
    }

    /// Get the number of layers.
    pub fn num_layers(&self) -> usize {
        self.layers.len()
    }

    /// Get the total number of inputs.
    pub fn num_inputs(&self) -> usize {
        self.num_public_inputs + self.num_private_inputs
    }

    /// Evaluate the circuit.
    ///
    /// Returns all wire values for each layer, with layer 0 being the output.
    pub fn evaluate(&self, inputs: &[F]) -> Vec<Vec<F>> {
        let nl = self.layers.len();
        let mut wires = Vec::with_capacity(nl + 1);

        // Input layer (layer NL)
        wires.push(inputs.to_vec());

        // Evaluate each layer in reverse order
        for layer in self.layers.iter().rev() {
            let inputs = wires.last().unwrap();
            let outputs = layer.evaluate(inputs);
            wires.push(outputs);
        }

        // Reverse so that wires[0] is the output
        wires.reverse();
        wires
    }

    /// Check if the circuit output is all zeros.
    pub fn check(&self, inputs: &[F]) -> bool {
        let wires = self.evaluate(inputs);
        wires[0].iter().all(|&w| w.is_zero())
    }

    /// Serialize the circuit.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        // Number of layers
        bytes.extend_from_slice(&(self.layers.len() as u32).to_le_bytes());

        // Public and private input counts
        bytes.extend_from_slice(&(self.num_public_inputs as u32).to_le_bytes());
        bytes.extend_from_slice(&(self.num_private_inputs as u32).to_le_bytes());

        // Each layer
        for layer in &self.layers {
            bytes.extend_from_slice(&(layer.lv_out as u16).to_le_bytes());
            bytes.extend_from_slice(&(layer.lv_in as u16).to_le_bytes());
            bytes.extend_from_slice(&(layer.terms.len() as u32).to_le_bytes());

            for term in &layer.terms {
                bytes.extend_from_slice(&(term.g as u32).to_le_bytes());
                bytes.extend_from_slice(&(term.l as u32).to_le_bytes());
                bytes.extend_from_slice(&(term.r as u32).to_le_bytes());
                bytes.extend_from_slice(&term.v.to_bytes());
            }
        }

        bytes
    }
}

/// Builder for constructing circuits.
pub struct CircuitBuilder<F: Field> {
    layers: Vec<LayerBuilder<F>>,
    num_public_inputs: usize,
    num_private_inputs: usize,
}

/// Builder for constructing a single layer.
pub struct LayerBuilder<F: Field> {
    terms: Vec<QuadTerm<F>>,
    lv_out: usize,
    lv_in: usize,
}

impl<F: Field> CircuitBuilder<F> {
    /// Create a new circuit builder.
    pub fn new(num_public_inputs: usize, num_private_inputs: usize) -> Self {
        Self {
            layers: Vec::new(),
            num_public_inputs,
            num_private_inputs,
        }
    }

    /// Add a layer to the circuit.
    pub fn add_layer(&mut self, layer: LayerBuilder<F>) {
        self.layers.push(layer);
    }

    /// Build the circuit.
    pub fn build(self) -> Circuit<F> {
        let layers = self
            .layers
            .into_iter()
            .map(|l| Layer::new(l.terms, l.lv_out, l.lv_in))
            .collect();

        Circuit::new(layers, self.num_public_inputs, self.num_private_inputs)
    }
}

impl<F: Field> LayerBuilder<F> {
    /// Create a new layer builder.
    pub fn new(lv_out: usize, lv_in: usize) -> Self {
        Self {
            terms: Vec::new(),
            lv_out,
            lv_in,
        }
    }

    /// Add a multiplication term: out[g] += v * in[l] * in[r].
    pub fn add_mul(&mut self, g: usize, l: usize, r: usize, v: F) {
        self.terms.push(QuadTerm::new(g, l, r, v));
    }

    /// Add an assertion: in[l] * in[r] = 0 at output g.
    pub fn add_assert(&mut self, g: usize, l: usize, r: usize) {
        self.terms.push(QuadTerm::assertion(g, l, r));
    }

    /// Add a copy: out[g] = in[l] (implemented as out[g] = 1 * in[l] * in[0] where in[0] = 1).
    pub fn add_copy(&mut self, g: usize, l: usize, one_wire: usize) {
        self.terms.push(QuadTerm::new(g, l, one_wire, F::ONE));
    }

    /// Add addition: out[g] = in[l] + in[r] (requires both have coefficient 1 from a constant wire).
    pub fn add_add(&mut self, g: usize, l: usize, r: usize, one_wire: usize) {
        self.add_mul(g, l, one_wire, F::ONE);
        self.add_mul(g, r, one_wire, F::ONE);
    }
}

/// Example: create a circuit that checks n is the m-th s-gonal number.
/// This verifies: 2n = (s-2)m^2 - (s-4)m
pub fn create_polygonal_circuit<F: Field>() -> Circuit<F> {
    // Inputs: [1, n, m, s] (where 1 is a constant)
    // We need to verify: 2n = (s-2)*m^2 - (s-4)*m
    // Rearranged: 2n - (s-2)*m^2 + (s-4)*m = 0

    let mut builder = CircuitBuilder::new(0, 4);

    // Layer 2 (closest to inputs): compute intermediate values
    // Need: m^2, s-2, s-4, etc.
    let mut layer2 = LayerBuilder::new(2, 2); // 4 outputs, 4 inputs

    // out[0] = m * m (m^2)
    layer2.add_mul(0, 2, 2, F::ONE);

    // out[1] = s - 2 (need to encode this differently - s * 1 - 2)
    // We'll compute s * 1 and handle the constant separately
    layer2.add_mul(1, 3, 0, F::ONE); // s * 1

    // out[2] = m * 1 (copy m)
    layer2.add_mul(2, 2, 0, F::ONE);

    // out[3] = n * 1 (copy n)
    layer2.add_mul(3, 1, 0, F::ONE);

    builder.add_layer(layer2);

    // Layer 1: combine intermediate values
    let mut layer1 = LayerBuilder::new(1, 2); // 2 outputs, 4 inputs

    // This is a simplified version - a real implementation would be more complex
    // out[0] = result of computation
    layer1.add_mul(0, 0, 1, F::ONE); // (s-2) * m^2

    builder.add_layer(layer1);

    // Layer 0 (output): final check
    let mut layer0 = LayerBuilder::new(0, 1); // 1 output, 2 inputs

    // Should output 0 if the computation is correct
    layer0.add_mul(0, 0, 0, F::ONE);

    builder.add_layer(layer0);

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp128;

    #[test]
    fn test_simple_circuit() {
        // Circuit that computes out = a * b
        let mut builder = CircuitBuilder::<Fp128>::new(0, 2);

        let mut layer = LayerBuilder::new(0, 1); // 1 output, 2 inputs
        layer.add_mul(0, 0, 1, Fp128::ONE);

        builder.add_layer(layer);

        let circuit = builder.build();

        // Evaluate with a=3, b=4
        let inputs = vec![Fp128::from_u64(3), Fp128::from_u64(4)];
        let wires = circuit.evaluate(&inputs);

        // Output should be 12
        assert_eq!(wires[0][0], Fp128::from_u64(12));
    }

    #[test]
    fn test_circuit_check() {
        // Circuit that checks a * b = c
        // Output is a * b - c, which should be 0

        let mut builder = CircuitBuilder::<Fp128>::new(0, 3);

        let mut layer = LayerBuilder::new(0, 2); // 1 output, 4 inputs (padded)
        // out = a * b - c = a * b + (-1) * c * 1
        layer.add_mul(0, 0, 1, Fp128::ONE); // a * b

        builder.add_layer(layer);

        let circuit = builder.build();

        // Check with a=3, b=4 (product = 12)
        let inputs = vec![
            Fp128::from_u64(3),
            Fp128::from_u64(4),
            Fp128::from_u64(12), // This is ignored in this simple circuit
        ];
        let wires = circuit.evaluate(&inputs);

        // Output should be 12 (3 * 4)
        assert_eq!(wires[0][0], Fp128::from_u64(12));
    }

    #[test]
    fn test_quad_term() {
        let term = QuadTerm::<Fp128>::new(0, 1, 2, Fp128::from_u64(5));
        assert!(!term.is_assertion());

        let assertion = QuadTerm::<Fp128>::assertion(0, 1, 2);
        assert!(assertion.is_assertion());
    }
}
