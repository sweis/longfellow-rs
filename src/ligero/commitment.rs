//! Ligero commitment structure.
//!
//! A Ligero commitment is the root of a Merkle tree formed from columns
//! of the tableau matrix.

use crate::merkle::MerkleDigest;

/// A Ligero commitment to a witness vector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LigeroCommitment {
    /// The Merkle root of the tableau columns.
    pub root: MerkleDigest,
}

impl LigeroCommitment {
    /// Create a new commitment from a Merkle root.
    pub fn new(root: MerkleDigest) -> Self {
        Self { root }
    }

    /// Serialize the commitment to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        self.root.to_vec()
    }

    /// Deserialize a commitment from bytes.
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 32 {
            return None;
        }
        let mut root = [0u8; 32];
        root.copy_from_slice(bytes);
        Some(Self { root })
    }
}
