//! Sumcheck protocol implementation.
//!
//! This module implements the sumcheck protocol used in the Longfellow ZK scheme
//! for verifiable computation on layered circuits.

mod eq;
mod layer;
mod proof;
mod verifier;

pub use eq::{bind, bind_eq, bindv};
pub use layer::{sumcheck_layer, LayerPad, PolyEvals, SumcheckLayerProof};
pub use proof::{sumcheck_circuit, CircuitLayer, CircuitPad, LayerWires, SumcheckProof};
pub use verifier::{
    verify_sumcheck_circuit, verify_sumcheck_layer, verify_sumcheck_rounds, LayerVerifyResult,
};
