//! Ligero commitment scheme implementation.
//!
//! This module implements the Ligero commitment scheme as described in the
//! Longfellow ZK protocol specification.

mod commitment;
mod params;
mod proof;
mod prover;
mod verifier;

pub use commitment::LigeroCommitment;
pub use params::LigeroParams;
pub use proof::LigeroProof;
pub use prover::{LigeroProver, LinearConstraintTerm, QuadraticConstraint};
pub use verifier::verify_ligero;
