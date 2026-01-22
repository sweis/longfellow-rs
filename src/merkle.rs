//! Merkle tree implementation for Ligero commitments.
//!
//! This module implements a Merkle tree with batch proof support as specified
//! in the Longfellow ZK protocol.

use sha2::{Digest, Sha256};

/// A 32-byte digest (SHA-256 output).
pub type MerkleDigest = [u8; 32];

/// A Merkle tree with n leaves.
///
/// The tree is stored in an array where:
/// - Index 0 is unused
/// - Index 1 is the root
/// - Indices 2..n are internal nodes
/// - Indices n..2n are leaves
#[derive(Clone, Debug)]
pub struct MerkleTree {
    /// The number of leaves.
    n: usize,
    /// Array of 2*n digests.
    nodes: Vec<MerkleDigest>,
}

impl MerkleTree {
    /// Create a new Merkle tree with the given number of leaves.
    pub fn new(n: usize) -> Self {
        assert!(n > 0, "Merkle tree must have at least one leaf");
        Self {
            n,
            nodes: vec![[0u8; 32]; 2 * n],
        }
    }

    /// Set a leaf at the given position (0-indexed).
    pub fn set_leaf(&mut self, pos: usize, leaf: MerkleDigest) {
        assert!(pos < self.n, "Leaf position out of bounds");
        self.nodes[pos + self.n] = leaf;
    }

    /// Build the tree and return the root.
    ///
    /// This must be called after all leaves are set.
    pub fn build(&mut self) -> MerkleDigest {
        // Build from leaves up to root
        for i in (1..self.n).rev() {
            let left = &self.nodes[2 * i];
            let right = &self.nodes[2 * i + 1];
            self.nodes[i] = hash_pair(left, right);
        }
        self.root()
    }

    /// Get the root of the tree.
    pub fn root(&self) -> MerkleDigest {
        self.nodes[1]
    }

    /// Get the number of leaves.
    pub fn num_leaves(&self) -> usize {
        self.n
    }

    /// Create a compressed batch proof for the given leaf indices.
    ///
    /// The indices should be 0-indexed positions of leaves to prove.
    pub fn compressed_proof(&self, requested_leaves: &[usize]) -> MerkleProof {
        let marked = self.mark_tree(requested_leaves);
        let mut proof_nodes = Vec::new();

        // Traverse from bottom to top (n-1 down to 1)
        for i in (1..self.n).rev() {
            if marked[i] {
                let child = 2 * i;
                // Determine which child is marked
                let sibling = if marked[child] { child + 1 } else { child };

                // If sibling is not marked, include it in the proof
                if !marked[sibling] {
                    proof_nodes.push(self.nodes[sibling]);
                }
            }
        }

        MerkleProof {
            nodes: proof_nodes,
            leaf_count: self.n,
        }
    }

    /// Mark the tree nodes that are on the path from leaves to root.
    fn mark_tree(&self, requested_leaves: &[usize]) -> Vec<bool> {
        let mut marked = vec![false; 2 * self.n];

        // Mark the requested leaves
        for &index in requested_leaves {
            marked[index + self.n] = true;
        }

        // Mark ancestors from bottom to top
        for i in (1..self.n).rev() {
            // Mark parent if any child is marked
            marked[i] = marked[2 * i] || marked[2 * i + 1];
        }

        marked
    }
}

/// A compressed Merkle proof for batch verification.
#[derive(Clone, Debug)]
pub struct MerkleProof {
    /// The nodes needed for verification.
    pub nodes: Vec<MerkleDigest>,
    /// The number of leaves in the tree.
    pub leaf_count: usize,
}

impl MerkleProof {
    /// Serialize the proof to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&(self.leaf_count as u64).to_le_bytes());
        bytes.extend_from_slice(&(self.nodes.len() as u64).to_le_bytes());
        for node in &self.nodes {
            bytes.extend_from_slice(node);
        }
        bytes
    }

    /// Deserialize a proof from bytes.
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() < 16 {
            return None;
        }

        let leaf_count = u64::from_le_bytes(bytes[0..8].try_into().ok()?) as usize;
        let node_count = u64::from_le_bytes(bytes[8..16].try_into().ok()?) as usize;

        if bytes.len() != 16 + node_count * 32 {
            return None;
        }

        let mut nodes = Vec::with_capacity(node_count);
        for i in 0..node_count {
            let start = 16 + i * 32;
            let mut node = [0u8; 32];
            node.copy_from_slice(&bytes[start..start + 32]);
            nodes.push(node);
        }

        Some(Self { nodes, leaf_count })
    }
}

/// Verify a compressed Merkle proof.
///
/// Verifies that the given leaf digests at the specified indices
/// belong to a tree with the given root.
pub fn verify_merkle(
    root: &MerkleDigest,
    n: usize,
    leaf_digests: &[MerkleDigest],
    indices: &[usize],
    proof: &MerkleProof,
) -> bool {
    if indices.is_empty() {
        return false;
    }

    if indices.len() != leaf_digests.len() {
        return false;
    }

    if n != proof.leaf_count {
        return false;
    }

    // Check for duplicate indices
    let mut sorted_indices = indices.to_vec();
    sorted_indices.sort();
    for i in 1..sorted_indices.len() {
        if sorted_indices[i] == sorted_indices[i - 1] {
            return false;
        }
    }

    // Build temporary arrays for reconstruction
    let mut tmp = vec![[0u8; 32]; 2 * n];
    let mut defined = vec![false; 2 * n];

    // Mark the tree
    let marked = mark_tree_verify(indices, n);

    // Read proof nodes
    let mut proof_index = 0;
    for i in (1..n).rev() {
        if marked[i] {
            let child = 2 * i;
            let sibling = if marked[child] { child + 1 } else { child };

            if !marked[sibling] {
                if proof_index >= proof.nodes.len() {
                    return false;
                }
                tmp[sibling] = proof.nodes[proof_index];
                defined[sibling] = true;
                proof_index += 1;
            }
        }
    }

    // Set the leaf values
    for (i, &index) in indices.iter().enumerate() {
        tmp[index + n] = leaf_digests[i];
        defined[index + n] = true;
    }

    // Recompute internal nodes from leaves to root
    for i in (1..n).rev() {
        if defined[2 * i] && defined[2 * i + 1] {
            tmp[i] = hash_pair(&tmp[2 * i], &tmp[2 * i + 1]);
            defined[i] = true;
        }
    }

    // Check root
    defined[1] && tmp[1] == *root
}

/// Mark tree nodes for verification.
fn mark_tree_verify(indices: &[usize], n: usize) -> Vec<bool> {
    let mut marked = vec![false; 2 * n];

    for &index in indices {
        marked[index + n] = true;
    }

    for i in (1..n).rev() {
        marked[i] = marked[2 * i] || marked[2 * i + 1];
    }

    marked
}

/// Hash two digests together.
pub fn hash_pair(left: &MerkleDigest, right: &MerkleDigest) -> MerkleDigest {
    let mut hasher = Sha256::new();
    hasher.update(left);
    hasher.update(right);
    hasher.finalize().into()
}

/// Hash arbitrary data to a digest.
pub fn hash_data(data: &[u8]) -> MerkleDigest {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merkle_basic() {
        let mut tree = MerkleTree::new(4);

        // Set leaves
        tree.set_leaf(0, hash_data(b"leaf0"));
        tree.set_leaf(1, hash_data(b"leaf1"));
        tree.set_leaf(2, hash_data(b"leaf2"));
        tree.set_leaf(3, hash_data(b"leaf3"));

        let root = tree.build();

        // Verify a proof for leaf 0
        let proof = tree.compressed_proof(&[0]);
        assert!(verify_merkle(
            &root,
            4,
            &[hash_data(b"leaf0")],
            &[0],
            &proof
        ));
    }

    #[test]
    fn test_merkle_batch_proof() {
        let mut tree = MerkleTree::new(8);

        for i in 0..8 {
            tree.set_leaf(i, hash_data(&[i as u8]));
        }

        let root = tree.build();

        // Batch proof for leaves 1, 3, 5
        let indices = vec![1, 3, 5];
        let leaves: Vec<_> = indices.iter().map(|&i| hash_data(&[i as u8])).collect();
        let proof = tree.compressed_proof(&indices);

        assert!(verify_merkle(&root, 8, &leaves, &indices, &proof));
    }

    #[test]
    fn test_merkle_invalid_proof() {
        let mut tree = MerkleTree::new(4);

        for i in 0..4 {
            tree.set_leaf(i, hash_data(&[i as u8]));
        }

        let root = tree.build();
        let proof = tree.compressed_proof(&[0]);

        // Try to verify with wrong leaf
        assert!(!verify_merkle(
            &root,
            4,
            &[hash_data(b"wrong")],
            &[0],
            &proof
        ));
    }

    #[test]
    fn test_merkle_spec_vector_1() {
        // From the spec test vectors
        let leaves = vec![
            hex::decode("4bf5122f344554c53bde2ebb8cd2b7e3d1600ad631c385a5d7cce23c7785459a")
                .unwrap(),
            hex::decode("dbc1b4c900ffe48d575b5da5c638040125f65db0fe3e24494b76ea986457d986")
                .unwrap(),
            hex::decode("084fed08b978af4d7d196a7446a86b58009e636b611db16211b65a9aadff29c5")
                .unwrap(),
            hex::decode("e52d9c508c502347344d8c07ad91cbd6068afc75ff6292f062a09ca381c89e71")
                .unwrap(),
            hex::decode("e77b9a9ae9e30b0dbdb6f510a264ef9de781501d7b6b92ae89eb059c5ab743db")
                .unwrap(),
        ];

        let expected_root =
            hex::decode("f22f4501ffd3bdffcecc9e4cd6828a4479aeedd6aa484eb7c1f808ccf71c6e76")
                .unwrap();

        let mut tree = MerkleTree::new(5);
        for (i, leaf) in leaves.iter().enumerate() {
            let mut digest = [0u8; 32];
            digest.copy_from_slice(leaf);
            tree.set_leaf(i, digest);
        }

        let root = tree.build();
        assert_eq!(root.to_vec(), expected_root);
    }
}
