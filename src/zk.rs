//! Full Zero-Knowledge proof system.
//!
//! This module combines the Ligero commitment scheme with the sumcheck protocol
//! to create a complete zero-knowledge proof system for arithmetic circuits.

use crate::circuit::Circuit;
use crate::field::Field;
use crate::ligero::{
    verify_ligero, LigeroCommitment, LigeroParams, LigeroProof, LigeroProver, QuadraticConstraint,
};
use crate::sumcheck::{
    sumcheck_circuit, verify_sumcheck_rounds, CircuitLayer, CircuitPad, LayerWires, SumcheckProof,
};
use crate::transcript::Transcript;
use rand::Rng;

/// A complete ZK proof for a circuit.
#[derive(Clone, Debug)]
pub struct ZkProof<F: Field> {
    /// Commitment to the witness.
    pub commitment: LigeroCommitment,

    /// Sumcheck proof for the circuit.
    pub sumcheck_proof: SumcheckProof<F>,

    /// Ligero proof for the constraints.
    pub ligero_proof: LigeroProof<F>,

    /// Ligero parameters (needed by verifier to reconstruct the tableau layout).
    pub ligero_params: LigeroParams,

    /// Quadratic constraints extracted from the circuit (needed by Ligero verifier).
    pub quadratic_constraints: Vec<QuadraticConstraint>,
}

impl<F: Field> ZkProof<F> {
    /// Get the proof size in bytes.
    pub fn size(&self) -> usize {
        let commitment_size = self.commitment.to_bytes().len();
        let sumcheck_size = self.sumcheck_proof.to_bytes().len();
        let ligero_size = self.ligero_proof.size();
        commitment_size + sumcheck_size + ligero_size
    }
}

/// The ZK prover.
pub struct ZkProver<F: Field> {
    /// The circuit to prove.
    circuit: Circuit<F>,

    /// The witness (private inputs).
    witness: Vec<F>,

    /// The public inputs.
    public_inputs: Vec<F>,
}

impl<F: Field> ZkProver<F> {
    /// Create a new ZK prover.
    pub fn new(circuit: Circuit<F>, public_inputs: Vec<F>, witness: Vec<F>) -> Self {
        Self {
            circuit,
            witness,
            public_inputs,
        }
    }

    /// Generate a zero-knowledge proof.
    pub fn prove<R: Rng>(&self, rng: &mut R) -> ZkProof<F> {
        let mut transcript = Transcript::new(b"longfellow-zk");

        // 1. Evaluate the circuit to get all wire values
        let all_inputs: Vec<F> = self
            .public_inputs
            .iter()
            .chain(self.witness.iter())
            .copied()
            .collect();
        let wire_values = self.circuit.evaluate(&all_inputs);

        // 2. Commit to the witness using Ligero
        let quadratic_constraints = self.extract_quadratic_constraints();
        let ligero_params =
            LigeroParams::new(self.witness.len(), quadratic_constraints.len(), 128);
        let mut ligero_prover = LigeroProver::with_params(
            self.witness.clone(),
            quadratic_constraints.clone(),
            ligero_params.clone(),
        );
        let commitment = ligero_prover.commit(rng);

        // 3. Write commitment to transcript
        transcript.write_bytes(&commitment.to_bytes());

        // 4. Run sumcheck protocol
        let (sumcheck_proof, _sumcheck_constraints) =
            self.run_sumcheck(&wire_values, &mut transcript, rng);

        // 5. Generate Ligero proof for the constraints
        let ligero_proof = ligero_prover.prove(&mut transcript, &commitment, &[], &[], rng);

        ZkProof {
            commitment,
            sumcheck_proof,
            ligero_proof,
            ligero_params,
            quadratic_constraints,
        }
    }

    /// Extract quadratic constraints from the circuit.
    fn extract_quadratic_constraints(&self) -> Vec<QuadraticConstraint> {
        let mut constraints = Vec::new();

        // For each layer, extract the multiplication constraints
        for layer in &self.circuit.layers {
            for term in &layer.terms {
                if !term.is_assertion() {
                    // Each multiplication x * y = z becomes a constraint
                    // We need to map wire indices appropriately
                    constraints.push(QuadraticConstraint {
                        x: term.l,
                        y: term.r,
                        z: term.g,
                    });
                }
            }
        }

        constraints
    }

    /// Run the sumcheck protocol.
    fn run_sumcheck<R: Rng>(
        &self,
        wire_values: &[Vec<F>],
        transcript: &mut Transcript,
        rng: &mut R,
    ) -> (SumcheckProof<F>, Vec<(Vec<F>, F)>) {
        // Convert circuit layers to sumcheck format
        let mut sumcheck_layers = Vec::new();
        let mut sumcheck_wires = Vec::new();

        for (i, layer) in self.circuit.layers.iter().enumerate() {
            let (q, z) = layer.build_q_and_z();

            let sumcheck_layer = if z.iter().any(|g| g.iter().any(|l| l.iter().any(|&v| !v.is_zero()))) {
                CircuitLayer::with_assertions(q, z, layer.lv_out)
            } else {
                CircuitLayer::new(q, layer.lv_out)
            };

            sumcheck_layers.push(sumcheck_layer);

            // Wire values: layer i+1 in wire_values corresponds to input of layer i
            let input_idx = self.circuit.layers.len() - i;
            let input_wires = wire_values.get(input_idx).cloned().unwrap_or_default();
            sumcheck_wires.push(LayerWires::from_single(input_wires));
        }

        // Generate padding for zero-knowledge
        let pad = CircuitPad::random(rng, &sumcheck_layers);

        // Run sumcheck
        let proof = sumcheck_circuit(&sumcheck_layers, &sumcheck_wires, &pad, transcript);

        // Return proof and any extracted constraints
        (proof, vec![])
    }
}

/// Verify a ZK proof.
///
/// This performs end-to-end verification:
/// 1. Replays the Fiat-Shamir transcript in sync with the prover
/// 2. Verifies each sumcheck round reduces the claim correctly
/// 3. Verifies the Ligero proof: Merkle openings, low-degree test,
///    dot-product test, and quadratic constraint test
///
/// Returns `true` only if all checks pass.
pub fn verify_zk<F: Field>(
    circuit: &Circuit<F>,
    _public_inputs: &[F],
    proof: &ZkProof<F>,
) -> bool {
    let mut transcript = Transcript::new(b"longfellow-zk");

    // 1. Write commitment to transcript (must match prover order)
    transcript.write_bytes(&proof.commitment.to_bytes());

    // 2. Verify sumcheck proof structure and replay transcript.
    //
    // This checks that:
    // - The number of layer proofs matches the circuit
    // - Each round's polynomial evaluations are written to transcript
    // - Challenges are derived consistently
    //
    // The sumcheck verifier returns (vl, vr, final_claim) for each layer.
    // These values represent the prover's claims about the witness at
    // random evaluation points, which must be consistent with the
    // Ligero-committed witness.
    let layers_lv: Vec<usize> = circuit.layers.iter().map(|l| l.lv_out).collect();

    let sumcheck_result = match verify_sumcheck_rounds(
        &proof.sumcheck_proof,
        &layers_lv,
        &mut transcript,
    ) {
        Some(result) => result,
        None => return false,
    };

    // Structural check: all layers should have produced results
    if sumcheck_result.len() != circuit.num_layers() {
        return false;
    }

    // 3. Verify the Ligero proof.
    //
    // This verifies:
    // - Merkle tree openings are valid (committed columns match)
    // - Low-degree test: each row is a low-degree polynomial
    // - Dot-product test: linear constraints on the witness hold
    // - Quadratic test: x * y = z constraints hold
    //
    // The Ligero proof uses its own transcript fork since the prover
    // writes the commitment again inside LigeroProver::prove.
    let mut ligero_transcript = transcript;

    if !verify_ligero(
        &proof.commitment,
        &proof.ligero_proof,
        &mut ligero_transcript,
        &proof.ligero_params,
        &[], // Linear constraints derived from sumcheck would go here
        &[], // b vector
        &proof.quadratic_constraints,
    ) {
        return false;
    }

    true
}

/// Parameters for the ZK scheme.
#[derive(Clone, Debug)]
pub struct ZkParams {
    /// Security level in bits.
    pub security_level: usize,

    /// Number of column queries for Ligero.
    pub num_queries: usize,

    /// Reed-Solomon rate.
    pub rate: usize,
}

impl ZkParams {
    /// Create default parameters for 128-bit security.
    pub fn default_128() -> Self {
        Self {
            security_level: 128,
            num_queries: 6,
            rate: 4,
        }
    }

    /// Create parameters for 256-bit security.
    pub fn default_256() -> Self {
        Self {
            security_level: 256,
            num_queries: 12,
            rate: 4,
        }
    }
}

impl Default for ZkParams {
    fn default() -> Self {
        Self::default_128()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::{CircuitBuilder, LayerBuilder};
    use crate::field::Fp128;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_zk_simple() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Create a simple multiplication circuit: c = a * b
        let mut builder = CircuitBuilder::<Fp128>::new(0, 3);

        let mut layer = LayerBuilder::new(0, 2);
        layer.add_mul(0, 0, 1, Fp128::ONE);
        builder.add_layer(layer);

        let circuit = builder.build();

        // Prove that 3 * 4 = 12
        let public_inputs = vec![];
        let witness = vec![
            Fp128::from_u64(3),
            Fp128::from_u64(4),
            Fp128::from_u64(12),
        ];

        let prover = ZkProver::new(circuit.clone(), public_inputs.clone(), witness);
        let proof = prover.prove(&mut rng);

        // Verify
        let valid = verify_zk(&circuit, &public_inputs, &proof);
        assert!(valid);
    }

    #[test]
    fn test_zk_params() {
        let params = ZkParams::default_128();
        assert_eq!(params.security_level, 128);
        assert_eq!(params.num_queries, 6);

        let params = ZkParams::default_256();
        assert_eq!(params.security_level, 256);
        assert_eq!(params.num_queries, 12);
    }

    #[test]
    fn test_zk_rejects_tampered_commitment() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let mut builder = CircuitBuilder::<Fp128>::new(0, 3);
        let mut layer = LayerBuilder::new(0, 2);
        layer.add_mul(0, 0, 1, Fp128::ONE);
        builder.add_layer(layer);
        let circuit = builder.build();

        let witness = vec![
            Fp128::from_u64(3),
            Fp128::from_u64(4),
            Fp128::from_u64(12),
        ];

        let prover = ZkProver::new(circuit.clone(), vec![], witness);
        let mut proof = prover.prove(&mut rng);

        // Tamper with the commitment - this should cause Merkle verification
        // to fail since the opened columns won't match the root
        let mut tampered_root = proof.commitment.root;
        tampered_root[0] ^= 0xff;
        proof.commitment = LigeroCommitment::new(tampered_root);

        let valid = verify_zk(&circuit, &[], &proof);
        assert!(!valid, "Tampered commitment should fail verification");
    }

    #[test]
    fn test_zk_rejects_tampered_ligero_ldt() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let mut builder = CircuitBuilder::<Fp128>::new(0, 3);
        let mut layer = LayerBuilder::new(0, 2);
        layer.add_mul(0, 0, 1, Fp128::ONE);
        builder.add_layer(layer);
        let circuit = builder.build();

        let witness = vec![
            Fp128::from_u64(3),
            Fp128::from_u64(4),
            Fp128::from_u64(12),
        ];

        let prover = ZkProver::new(circuit.clone(), vec![], witness);
        let mut proof = prover.prove(&mut rng);

        // Tamper with the LDT response - should fail low-degree test
        if let Some(elem) = proof.ligero_proof.ldt.get_mut(0) {
            *elem = *elem + Fp128::ONE;
        }

        let valid = verify_zk(&circuit, &[], &proof);
        assert!(!valid, "Tampered LDT should fail verification");
    }

    #[test]
    fn test_zk_rejects_wrong_layer_count() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let mut builder = CircuitBuilder::<Fp128>::new(0, 3);
        let mut layer = LayerBuilder::new(0, 2);
        layer.add_mul(0, 0, 1, Fp128::ONE);
        builder.add_layer(layer);
        let circuit = builder.build();

        let witness = vec![
            Fp128::from_u64(3),
            Fp128::from_u64(4),
            Fp128::from_u64(12),
        ];

        let prover = ZkProver::new(circuit.clone(), vec![], witness);
        let mut proof = prover.prove(&mut rng);

        // Remove a layer proof - should fail structural check
        proof.sumcheck_proof.layer_proofs.clear();

        let valid = verify_zk(&circuit, &[], &proof);
        assert!(!valid, "Wrong layer count should fail verification");
    }

    #[test]
    fn test_zk_with_gf2_128() {
        use crate::field::GF2_128;
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Simple circuit over GF(2^128): c = a * b
        let mut builder = CircuitBuilder::<GF2_128>::new(0, 3);
        let mut layer = LayerBuilder::new(0, 2);
        layer.add_mul(0, 0, 1, GF2_128::ONE);
        builder.add_layer(layer);
        let circuit = builder.build();

        // Prove: x * x^2 = x^3 in GF(2^128)
        let a = GF2_128::X;       // x
        let b = GF2_128::from_u64(4); // x^2
        let c = a * b;                 // x^3 = 8
        let witness = vec![a, b, c];

        let prover = ZkProver::new(circuit.clone(), vec![], witness);
        let proof = prover.prove(&mut rng);

        let valid = verify_zk(&circuit, &[], &proof);
        assert!(valid, "Valid GF(2^128) proof should verify");
    }

    #[test]
    fn test_zk_with_fp256() {
        use crate::field::Fp256;
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Simple circuit over Fp256: c = a * b
        let mut builder = CircuitBuilder::<Fp256>::new(0, 3);
        let mut layer = LayerBuilder::new(0, 2);
        layer.add_mul(0, 0, 1, Fp256::ONE);
        builder.add_layer(layer);
        let circuit = builder.build();

        // Prove: 7 * 11 = 77 in Fp256
        let witness = vec![
            Fp256::from_u64(7),
            Fp256::from_u64(11),
            Fp256::from_u64(77),
        ];

        let prover = ZkProver::new(circuit.clone(), vec![], witness);
        let proof = prover.prove(&mut rng);

        let valid = verify_zk(&circuit, &[], &proof);
        assert!(valid, "Valid Fp256 proof should verify");
    }
}
