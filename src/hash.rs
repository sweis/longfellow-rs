//! Hash function abstraction for the Longfellow ZK scheme.
//!
//! This module provides a trait for hash functions that can be used throughout
//! the library, with implementations for SHA-256 and optionally BLAKE3.
//!
//! # Features
//!
//! - `default` - Uses SHA-256 (required by the spec)
//! - `blake3_hash` - Enables BLAKE3 support for faster hashing

use sha2::{Digest, Sha256};

/// A 32-byte hash digest.
pub type HashDigest = [u8; 32];

/// Trait for hash functions used in the Longfellow ZK scheme.
///
/// All hash functions must produce 32-byte digests to be compatible
/// with the protocol.
pub trait HashFunction: Clone + Default {
    /// The name of the hash function (for debugging/logging).
    const NAME: &'static str;

    /// Create a new hasher instance.
    fn new() -> Self;

    /// Update the hasher with data.
    fn update(&mut self, data: &[u8]);

    /// Finalize and return the hash digest.
    fn finalize(self) -> HashDigest;

    /// Convenience method to hash data in one call.
    fn hash(data: &[u8]) -> HashDigest {
        let mut h = Self::new();
        h.update(data);
        h.finalize()
    }

    /// Hash two digests together (for Merkle trees).
    fn hash_pair(left: &HashDigest, right: &HashDigest) -> HashDigest {
        let mut h = Self::new();
        h.update(left);
        h.update(right);
        h.finalize()
    }
}

/// SHA-256 hash function implementation.
#[derive(Clone, Default)]
pub struct Sha256Hash {
    hasher: Sha256,
}

impl HashFunction for Sha256Hash {
    const NAME: &'static str = "SHA-256";

    fn new() -> Self {
        Self {
            hasher: Sha256::new(),
        }
    }

    fn update(&mut self, data: &[u8]) {
        self.hasher.update(data);
    }

    fn finalize(self) -> HashDigest {
        self.hasher.finalize().into()
    }
}

/// BLAKE3 hash function implementation.
///
/// BLAKE3 is significantly faster than SHA-256 while maintaining
/// strong security properties. It is available when the `blake3_hash`
/// feature is enabled.
#[cfg(feature = "blake3_hash")]
#[derive(Clone, Default)]
pub struct Blake3Hash {
    hasher: blake3::Hasher,
}

#[cfg(feature = "blake3_hash")]
impl HashFunction for Blake3Hash {
    const NAME: &'static str = "BLAKE3";

    fn new() -> Self {
        Self {
            hasher: blake3::Hasher::new(),
        }
    }

    fn update(&mut self, data: &[u8]) {
        self.hasher.update(data);
    }

    fn finalize(self) -> HashDigest {
        self.hasher.finalize().into()
    }
}

/// The default hash function type.
///
/// Uses SHA-256 by default (as required by the spec), but can be overridden
/// to BLAKE3 by enabling the `blake3_hash` feature.
#[cfg(not(feature = "blake3_hash"))]
pub type DefaultHash = Sha256Hash;

/// The default hash function type (BLAKE3 when feature is enabled).
///
/// When the `blake3_hash` feature is enabled, BLAKE3 is used as the default
/// hash function for improved performance.
#[cfg(feature = "blake3_hash")]
pub type DefaultHash = Blake3Hash;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sha256_hash() {
        let digest = Sha256Hash::hash(b"hello world");
        assert_eq!(digest.len(), 32);

        // Known SHA-256 hash of "hello world"
        let expected = hex::decode("b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9")
            .unwrap();
        assert_eq!(digest.to_vec(), expected);
    }

    #[test]
    fn test_sha256_hash_pair() {
        let left = Sha256Hash::hash(b"left");
        let right = Sha256Hash::hash(b"right");
        let combined = Sha256Hash::hash_pair(&left, &right);

        // Should be deterministic
        let combined2 = Sha256Hash::hash_pair(&left, &right);
        assert_eq!(combined, combined2);

        // Order matters
        let reversed = Sha256Hash::hash_pair(&right, &left);
        assert_ne!(combined, reversed);
    }

    #[cfg(feature = "blake3_hash")]
    #[test]
    fn test_blake3_hash() {
        let digest = Blake3Hash::hash(b"hello world");
        assert_eq!(digest.len(), 32);

        // Known BLAKE3 hash of "hello world"
        let expected = hex::decode("d74981efa70a0c880b8d8c1985d075dbcbf679b99a5f9914e5aaf96b831a9e24")
            .unwrap();
        assert_eq!(digest.to_vec(), expected);
    }

    #[cfg(feature = "blake3_hash")]
    #[test]
    fn test_blake3_hash_pair() {
        let left = Blake3Hash::hash(b"left");
        let right = Blake3Hash::hash(b"right");
        let combined = Blake3Hash::hash_pair(&left, &right);

        // Should be deterministic
        let combined2 = Blake3Hash::hash_pair(&left, &right);
        assert_eq!(combined, combined2);
    }

    #[test]
    fn test_default_hash() {
        // DefaultHash should work regardless of features
        let digest = DefaultHash::hash(b"test");
        assert_eq!(digest.len(), 32);
    }
}
