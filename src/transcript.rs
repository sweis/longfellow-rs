//! Fiat-Shamir transcript for non-interactive proofs.
//!
//! This module implements the Fiat-Shamir heuristic to transform an interactive
//! proof protocol into a non-interactive one using a cryptographic hash function.
//!
//! # Hash Functions
//!
//! By default, the transcript uses SHA-256 as required by the spec.
//! When the `blake3_hash` feature is enabled, BLAKE3 is used instead
//! for improved performance.

use crate::field::Field;
use crate::hash::{DefaultHash, HashDigest, HashFunction};
use std::marker::PhantomData;

/// A Fiat-Shamir transcript generic over the hash function.
///
/// The transcript accumulates all messages exchanged between prover and verifier,
/// and uses a cryptographic hash function to generate pseudorandom challenges.
#[derive(Clone)]
pub struct TranscriptGeneric<H: HashFunction> {
    /// The current hash state (accumulated data).
    state: Vec<u8>,
    /// Buffer for generating multiple challenges from a single hash.
    challenge_buffer: Vec<u8>,
    /// Current position in the challenge buffer.
    buffer_pos: usize,
    /// Phantom data for the hash function type.
    _hash: PhantomData<H>,
}

impl<H: HashFunction> TranscriptGeneric<H> {
    /// Create a new transcript with a domain separator.
    ///
    /// The domain separator ensures that transcripts for different protocols
    /// are distinct even if they have the same messages.
    pub fn new(domain_separator: &[u8]) -> Self {
        let mut state = Vec::new();
        // Write the length of the domain separator followed by its contents
        state.extend_from_slice(&(domain_separator.len() as u32).to_le_bytes());
        state.extend_from_slice(domain_separator);

        Self {
            state,
            challenge_buffer: Vec::new(),
            buffer_pos: 0,
            _hash: PhantomData,
        }
    }

    /// Write raw bytes to the transcript.
    pub fn write_bytes(&mut self, data: &[u8]) {
        // Invalidate the challenge buffer when new data is written
        self.challenge_buffer.clear();
        self.buffer_pos = 0;

        // Write length-prefixed data
        self.state.extend_from_slice(&(data.len() as u64).to_le_bytes());
        self.state.extend_from_slice(data);
    }

    /// Write a field element to the transcript.
    pub fn write_field_element<F: Field>(&mut self, element: &F) {
        let bytes = element.to_bytes();
        self.write_bytes(&bytes);
    }

    /// Write an array of field elements to the transcript.
    pub fn write_field_elements<F: Field>(&mut self, elements: &[F]) {
        // Write the array length first
        self.challenge_buffer.clear();
        self.buffer_pos = 0;

        self.state.extend_from_slice(&(elements.len() as u64).to_le_bytes());
        for elem in elements {
            let bytes = elem.to_bytes();
            self.state.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
            self.state.extend_from_slice(&bytes);
        }
    }

    /// Generate a single field element challenge.
    pub fn generate_challenge<F: Field>(&mut self) -> F {
        self.generate_challenges(1)[0]
    }

    /// Generate multiple field element challenges.
    pub fn generate_challenges<F: Field>(&mut self, count: usize) -> Vec<F> {
        let mut challenges = Vec::with_capacity(count);

        // Determine the byte size needed for field elements
        let elem_size = F::ZERO.to_bytes().len();

        for _i in 0..count {
            // If we need more random bytes, extend the buffer
            while self.buffer_pos + elem_size > self.challenge_buffer.len() {
                self.extend_challenge_buffer();
            }

            // Extract bytes for this challenge
            let bytes = &self.challenge_buffer[self.buffer_pos..self.buffer_pos + elem_size];
            let challenge = F::from_bytes(bytes);
            challenges.push(challenge);

            self.buffer_pos += elem_size;
        }

        challenges
    }

    /// Extend the challenge buffer by hashing the current state.
    fn extend_challenge_buffer(&mut self) {
        // Add a counter to ensure different extensions produce different output
        let counter = (self.challenge_buffer.len() / 32) as u64;

        let mut hasher = H::new();
        hasher.update(&self.state);
        hasher.update(&counter.to_le_bytes());

        let digest: HashDigest = hasher.finalize();
        self.challenge_buffer.extend_from_slice(&digest);
    }

    /// Generate random bytes for sampling without replacement.
    fn generate_bytes(&mut self, count: usize) -> Vec<u8> {
        while self.buffer_pos + count > self.challenge_buffer.len() {
            self.extend_challenge_buffer();
        }

        let bytes = self.challenge_buffer[self.buffer_pos..self.buffer_pos + count].to_vec();
        self.buffer_pos += count;
        bytes
    }

    /// Generate `count` unique indices in the range [0, n) without replacement.
    ///
    /// This uses the Fisher-Yates shuffle algorithm with rejection sampling.
    pub fn generate_indices_without_replacement(&mut self, n: usize, count: usize) -> Vec<usize> {
        assert!(count <= n, "Cannot sample more elements than available");

        let mut result = Vec::with_capacity(count);
        let mut available: Vec<usize> = (0..n).collect();

        for i in 0..count {
            // Generate a random index in [0, n - i)
            let remaining = n - i;
            let idx = self.generate_random_index(remaining);

            // Swap and collect
            let chosen = available[idx];
            available[idx] = available[remaining - 1];
            result.push(chosen);
        }

        result
    }

    /// Generate a random index in [0, n).
    fn generate_random_index(&mut self, n: usize) -> usize {
        // Use rejection sampling for uniform distribution
        let bits_needed = (usize::BITS - n.leading_zeros()) as usize;
        let bytes_needed = (bits_needed + 7) / 8;
        let mask = (1usize << bits_needed) - 1;

        loop {
            let bytes = self.generate_bytes(bytes_needed);
            let mut value = 0usize;
            for (i, &b) in bytes.iter().enumerate() {
                value |= (b as usize) << (i * 8);
            }
            value &= mask;

            if value < n {
                return value;
            }
        }
    }

    /// Fork the transcript to create an independent copy.
    pub fn fork(&self) -> Self {
        self.clone()
    }
}

/// Type alias for the default Transcript (uses DefaultHash).
pub type Transcript = TranscriptGeneric<DefaultHash>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp128;
    use crate::hash::Sha256Hash;

    #[test]
    fn test_transcript_determinism() {
        let mut t1 = Transcript::new(b"test");
        let mut t2 = Transcript::new(b"test");

        t1.write_bytes(&[1, 2, 3]);
        t2.write_bytes(&[1, 2, 3]);

        let c1: Fp128 = t1.generate_challenge();
        let c2: Fp128 = t2.generate_challenge();

        assert_eq!(c1, c2);
    }

    #[test]
    fn test_different_domain_separators() {
        let mut t1 = Transcript::new(b"test1");
        let mut t2 = Transcript::new(b"test2");

        let c1: Fp128 = t1.generate_challenge();
        let c2: Fp128 = t2.generate_challenge();

        assert_ne!(c1, c2);
    }

    #[test]
    fn test_indices_without_replacement() {
        let mut t = Transcript::new(b"test");
        let indices = t.generate_indices_without_replacement(100, 10);

        assert_eq!(indices.len(), 10);

        // Check uniqueness
        let mut sorted = indices.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), 10);

        // Check range
        for &idx in &indices {
            assert!(idx < 100);
        }
    }

    #[test]
    fn test_write_field_elements() {
        let mut t = Transcript::new(b"test");

        let elements = vec![Fp128::from_u64(1), Fp128::from_u64(2), Fp128::from_u64(3)];
        t.write_field_elements(&elements);

        let challenge: Fp128 = t.generate_challenge();
        assert!(!challenge.is_zero());
    }

    #[test]
    fn test_sha256_transcript_explicit() {
        // Test using SHA-256 explicitly
        let mut t1 = TranscriptGeneric::<Sha256Hash>::new(b"test");
        let mut t2 = TranscriptGeneric::<Sha256Hash>::new(b"test");

        t1.write_bytes(&[1, 2, 3]);
        t2.write_bytes(&[1, 2, 3]);

        let c1: Fp128 = t1.generate_challenge();
        let c2: Fp128 = t2.generate_challenge();

        assert_eq!(c1, c2);
    }

    #[cfg(feature = "blake3_hash")]
    #[test]
    fn test_blake3_transcript() {
        use crate::hash::Blake3Hash;

        let mut t1 = TranscriptGeneric::<Blake3Hash>::new(b"test");
        let mut t2 = TranscriptGeneric::<Blake3Hash>::new(b"test");

        t1.write_bytes(&[1, 2, 3]);
        t2.write_bytes(&[1, 2, 3]);

        let c1: Fp128 = t1.generate_challenge();
        let c2: Fp128 = t2.generate_challenge();

        assert_eq!(c1, c2);
    }

    #[cfg(feature = "blake3_hash")]
    #[test]
    fn test_different_hash_different_output() {
        use crate::hash::Blake3Hash;

        let mut t_sha = TranscriptGeneric::<Sha256Hash>::new(b"test");
        let mut t_blake = TranscriptGeneric::<Blake3Hash>::new(b"test");

        t_sha.write_bytes(&[1, 2, 3]);
        t_blake.write_bytes(&[1, 2, 3]);

        let c_sha: Fp128 = t_sha.generate_challenge();
        let c_blake: Fp128 = t_blake.generate_challenge();

        // Different hash functions should produce different challenges
        assert_ne!(c_sha, c_blake);
    }
}
