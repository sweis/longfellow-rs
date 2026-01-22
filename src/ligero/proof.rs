//! Ligero proof structure.
//!
//! This module defines the structure of a Ligero proof, which includes
//! the responses to the low-degree test, linear test, quadratic test,
//! and the Merkle proof for the opened columns.

use crate::field::Field;
use crate::merkle::MerkleProof;

/// A Ligero proof for a committed witness.
#[derive(Clone, Debug)]
pub struct LigeroProof<F: Field> {
    /// Low-degree test response: first BLOCK elements of the combined row.
    pub ldt: Vec<F>,

    /// Dot product test response.
    pub dot: Vec<F>,

    /// Quadratic proof: full DBLOCK evaluations.
    pub qpr: Vec<F>,

    /// Opened columns at the challenged indices.
    /// Each column contains one element per row of the tableau.
    pub columns: Vec<Vec<F>>,

    /// Merkle proof for the opened columns.
    pub merkle_proof: MerkleProof,
}

impl<F: Field> LigeroProof<F> {
    /// Create a new Ligero proof.
    pub fn new(
        ldt: Vec<F>,
        dot: Vec<F>,
        qpr: Vec<F>,
        columns: Vec<Vec<F>>,
        merkle_proof: MerkleProof,
    ) -> Self {
        Self {
            ldt,
            dot,
            qpr,
            columns,
            merkle_proof,
        }
    }

    /// Serialize the proof to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();

        // Write LDT response
        bytes.extend_from_slice(&(self.ldt.len() as u32).to_le_bytes());
        for elem in &self.ldt {
            bytes.extend_from_slice(&elem.to_bytes());
        }

        // Write dot response
        bytes.extend_from_slice(&(self.dot.len() as u32).to_le_bytes());
        for elem in &self.dot {
            bytes.extend_from_slice(&elem.to_bytes());
        }

        // Write quadratic proof
        bytes.extend_from_slice(&(self.qpr.len() as u32).to_le_bytes());
        for elem in &self.qpr {
            bytes.extend_from_slice(&elem.to_bytes());
        }

        // Write columns
        bytes.extend_from_slice(&(self.columns.len() as u32).to_le_bytes());
        if !self.columns.is_empty() {
            bytes.extend_from_slice(&(self.columns[0].len() as u32).to_le_bytes());
            for col in &self.columns {
                for elem in col {
                    bytes.extend_from_slice(&elem.to_bytes());
                }
            }
        }

        // Write Merkle proof
        let merkle_bytes = self.merkle_proof.to_bytes();
        bytes.extend_from_slice(&(merkle_bytes.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&merkle_bytes);

        bytes
    }

    /// Get the proof size in bytes.
    pub fn size(&self) -> usize {
        self.to_bytes().len()
    }
}
