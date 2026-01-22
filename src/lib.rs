//! Longfellow ZK - A Rust implementation of the Longfellow Zero-Knowledge proof system.
//!
//! This library implements the Longfellow ZK scheme as described in the IETF draft
//! [draft-google-cfrg-libzk](https://datatracker.ietf.org/doc/draft-google-cfrg-libzk/).
//!
//! # Overview
//!
//! The Longfellow ZK scheme is a succinct non-interactive zero-knowledge argument system
//! that combines:
//!
//! - **Ligero**: A commitment scheme that supports efficient ZK arguments for linear
//!   and quadratic constraints using Reed-Solomon codes and Merkle trees.
//!
//! - **Sumcheck**: A protocol for verifiable computation that enables efficient
//!   verification of polynomial evaluations over layered arithmetic circuits.
//!
//! The scheme requires only a collision-resistant hash function (SHA-256) and no
//! trusted setup or common reference string.
//!
//! # Example
//!
//! ```rust
//! use longfellow_zk::{
//!     circuit::{CircuitBuilder, LayerBuilder},
//!     field::Fp128,
//!     zk::{ZkProver, verify_zk},
//! };
//! use rand::SeedableRng;
//! use rand_chacha::ChaCha20Rng;
//!
//! // Create a circuit that checks a * b = c
//! let mut builder = CircuitBuilder::<Fp128>::new(0, 3);
//! let mut layer = LayerBuilder::new(0, 2);
//! layer.add_mul(0, 0, 1, Fp128::ONE);
//! builder.add_layer(layer);
//! let circuit = builder.build();
//!
//! // Create a proof
//! let mut rng = ChaCha20Rng::seed_from_u64(42);
//! let witness = vec![Fp128::from_u64(3), Fp128::from_u64(4), Fp128::from_u64(12)];
//! let prover = ZkProver::new(circuit.clone(), vec![], witness);
//! let proof = prover.prove(&mut rng);
//!
//! // Verify the proof
//! assert!(verify_zk(&circuit, &[], &proof));
//! ```
//!
//! # Modules
//!
//! - [`field`]: Finite field arithmetic over Fp128 = 2^128 - 2^108 + 1
//! - [`transcript`]: Fiat-Shamir transcript for non-interactive proofs
//! - [`merkle`]: Merkle tree implementation for commitments
//! - [`polynomial`]: Polynomial operations and Reed-Solomon encoding
//! - [`ligero`]: Ligero commitment scheme
//! - [`sumcheck`]: Sumcheck protocol for verifiable computation
//! - [`circuit`]: Layered arithmetic circuit representation
//! - [`zk`]: Full zero-knowledge prover and verifier
//!
//! # Security
//!
//! This implementation targets 128-bit security by default. The security relies on:
//!
//! - Collision resistance of SHA-256
//! - Soundness of the Ligero commitment scheme
//! - Soundness of the sumcheck protocol
//!
//! # References
//!
//! - [IETF Draft: The Longfellow Zero-knowledge Scheme](https://datatracker.ietf.org/doc/draft-google-cfrg-libzk/)
//! - [ePrint: Anonymous credentials from ECDSA](https://eprint.iacr.org/2024/2010)
//! - [Ligero: Lightweight Sublinear Arguments Without a Trusted Setup](https://eprint.iacr.org/2022/1608)

#![warn(missing_docs)]
#![warn(rust_2018_idioms)]

pub mod circuit;
pub mod field;
pub mod hash;
pub mod ligero;
pub mod merkle;
pub mod polynomial;
pub mod sumcheck;
pub mod transcript;
pub mod zk;

// Re-export commonly used types
pub use circuit::{Circuit, CircuitBuilder, Layer, LayerBuilder, QuadTerm};
pub use field::{batch_invert, Field, Fp128};
pub use hash::{DefaultHash, HashDigest, HashFunction, Sha256Hash};
#[cfg(feature = "blake3_hash")]
pub use hash::Blake3Hash;
pub use ligero::{LigeroCommitment, LigeroParams, LigeroProof, LigeroProver};
pub use merkle::{MerkleDigest, MerkleProof, MerkleTree, MerkleTreeGeneric};
pub use sumcheck::SumcheckProof;
pub use transcript::Transcript;
pub use zk::{ZkParams, ZkProof, ZkProver, verify_zk};

/// The crate version.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Error types for the library.
pub mod error {
    use thiserror::Error;

    /// Errors that can occur during proof generation or verification.
    #[derive(Error, Debug)]
    pub enum ZkError {
        /// Invalid witness provided.
        #[error("invalid witness: {0}")]
        InvalidWitness(String),

        /// Circuit evaluation failed.
        #[error("circuit evaluation failed: {0}")]
        CircuitError(String),

        /// Proof verification failed.
        #[error("proof verification failed: {0}")]
        VerificationError(String),

        /// Serialization/deserialization error.
        #[error("serialization error: {0}")]
        SerializationError(String),

        /// Invalid parameters.
        #[error("invalid parameters: {0}")]
        InvalidParams(String),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_version() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn test_basic_workflow() {
        let mut rng = ChaCha20Rng::seed_from_u64(12345);

        // Create a simple circuit
        let mut builder = CircuitBuilder::<Fp128>::new(0, 2);
        let mut layer = LayerBuilder::new(0, 1);
        layer.add_mul(0, 0, 1, Fp128::ONE);
        builder.add_layer(layer);
        let circuit = builder.build();

        // Create a proof for 5 * 7 = 35
        let witness = vec![Fp128::from_u64(5), Fp128::from_u64(7)];
        let prover = ZkProver::new(circuit.clone(), vec![], witness);
        let proof = prover.prove(&mut rng);

        // Verify
        assert!(verify_zk(&circuit, &[], &proof));
    }

    #[test]
    fn test_field_operations() {
        let a = Fp128::from_u64(123);
        let b = Fp128::from_u64(456);

        // Test basic operations
        let sum = a + b;
        let _diff = a - b;
        let prod = a * b;

        assert_eq!(sum, Fp128::from_u64(579));
        assert_eq!(prod, Fp128::from_u64(56088));

        // Test inverse
        let a_inv = a.invert().unwrap();
        assert_eq!(a * a_inv, Fp128::ONE);
    }

    #[test]
    fn test_transcript_determinism() {
        let mut t1 = Transcript::new(b"test");
        let mut t2 = Transcript::new(b"test");

        t1.write_bytes(&[1, 2, 3, 4]);
        t2.write_bytes(&[1, 2, 3, 4]);

        let c1: Fp128 = t1.generate_challenge();
        let c2: Fp128 = t2.generate_challenge();

        assert_eq!(c1, c2);
    }
}
