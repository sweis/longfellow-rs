//! Finite field arithmetic for the Longfellow ZK scheme.
//!
//! This module provides implementations of the finite fields used in the
//! Longfellow ZK protocol as specified in the IETF draft-google-cfrg-libzk:
//!
//! | Field       | FieldID | Description                                       |
//! |-------------|---------|---------------------------------------------------|
//! | Fp256       | 0x01    | NIST P-256 base field (prime = 2^256-2^224+2^192+2^96-1) |
//! | GF2_128     | 0x04    | Binary extension field GF(2)[x]/(x^128+x^7+x^2+x+1) |
//! | Fp128       | 0x06    | Prime field with p = 2^128 - 2^108 + 1            |
//!
//! # Field Trait
//!
//! All fields implement the [`Field`] trait, making the proof system
//! components (Ligero, sumcheck) generic over the choice of field.
//!
//! # The P2 Constant
//!
//! The sumcheck protocol represents degree-2 polynomials by their evaluations
//! at three points P0=0, P1=1, and P2. The choice of P2 depends on the field:
//!
//! - For prime fields (characteristic > 2): P2 = 2
//! - For GF(2^128): P2 = X (the polynomial variable, represented as 0b10)

mod fp128;
mod fp256;
mod gf2_128;

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub use fp128::Fp128;
pub use fp256::Fp256;
pub use gf2_128::GF2_128;

/// Field identifier per the IETF draft-google-cfrg-libzk specification.
///
/// These IDs are used in serialization to identify which field a proof
/// was constructed over.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum FieldId {
    /// NIST P-256 base field (prime = 2^256 - 2^224 + 2^192 + 2^96 - 1).
    P256 = 0x01,
    /// NIST P-384 base field.
    P384 = 0x02,
    /// NIST P-521 base field.
    P521 = 0x03,
    /// Binary extension field GF(2^128) = GF(2)[x]/(x^128 + x^7 + x^2 + x + 1).
    Gf2_128 = 0x04,
    /// Binary extension field GF(2^16), a subfield of GF(2^128).
    Gf2_16 = 0x05,
    /// Prime field with modulus 2^128 - 2^108 + 1.
    Fp128 = 0x06,
    /// Prime field with modulus 2^64 - 59.
    Fp64m59 = 0x07,
    /// Prime field with modulus 2^64 - 2^32 + 1 (Goldilocks).
    Goldilocks = 0x08,
    /// Quadratic extension of Fp64m59.
    Fp64m59Sq = 0x09,
    /// secp256k1 base field.
    Secp256 = 0x0a,
}

impl FieldId {
    /// Get the FieldId as a byte for serialization.
    pub fn as_byte(self) -> u8 {
        self as u8
    }

    /// Parse a FieldId from a byte.
    pub fn from_byte(b: u8) -> Option<Self> {
        match b {
            0x01 => Some(Self::P256),
            0x02 => Some(Self::P384),
            0x03 => Some(Self::P521),
            0x04 => Some(Self::Gf2_128),
            0x05 => Some(Self::Gf2_16),
            0x06 => Some(Self::Fp128),
            0x07 => Some(Self::Fp64m59),
            0x08 => Some(Self::Goldilocks),
            0x09 => Some(Self::Fp64m59Sq),
            0x0a => Some(Self::Secp256),
            _ => None,
        }
    }
}

/// Trait for field elements used in the Longfellow ZK scheme.
pub trait Field:
    Sized
    + Clone
    + Copy
    + Default
    + PartialEq
    + Eq
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Neg<Output = Self>
    + std::iter::Sum
    + std::iter::Product
    + std::fmt::Debug
{
    /// The additive identity.
    const ZERO: Self;

    /// The multiplicative identity.
    const ONE: Self;

    /// The evaluation point P2 for polynomial representation in sumcheck.
    ///
    /// - For prime fields with characteristic > 2, this is 2.
    /// - For GF(2^128), this is X (the polynomial variable).
    const P2: Self;

    /// The field ID per the IETF spec.
    const FIELD_ID: FieldId;

    /// The number of bytes in the serialized representation.
    const BYTE_LEN: usize;

    /// The characteristic of the field (2 for binary fields, else the prime).
    /// Returns 0 for large primes that don't fit in u64.
    const CHARACTERISTIC: u64;

    /// Create from a u64 value.
    fn from_u64(value: u64) -> Self;

    /// Check if zero.
    fn is_zero(&self) -> bool;

    /// Compute multiplicative inverse.
    fn invert(&self) -> Option<Self>;

    /// Square the element.
    fn square(&self) -> Self;

    /// Generate a random element.
    fn random<R: rand::Rng>(rng: &mut R) -> Self;

    /// Convert to bytes (little-endian).
    fn to_bytes(&self) -> Vec<u8>;

    /// Create from bytes (little-endian).
    fn from_bytes(bytes: &[u8]) -> Self;
}

/// Batch inversion using Montgomery's trick.
///
/// Given a slice of field elements [a_0, a_1, ..., a_{n-1}], computes
/// [a_0^{-1}, a_1^{-1}, ..., a_{n-1}^{-1}] using only one field inversion
/// and 3(n-1) multiplications.
///
/// Returns None if any element is zero.
pub fn batch_invert<F: Field>(elements: &[F]) -> Option<Vec<F>> {
    let n = elements.len();
    if n == 0 {
        return Some(vec![]);
    }

    // Check for zeros
    for e in elements {
        if e.is_zero() {
            return None;
        }
    }

    if n == 1 {
        return Some(vec![elements[0].invert()?]);
    }

    // Compute prefix products: prefix[i] = a_0 * a_1 * ... * a_i
    let mut prefix = Vec::with_capacity(n);
    prefix.push(elements[0]);
    for i in 1..n {
        prefix.push(prefix[i - 1] * elements[i]);
    }

    // Invert the total product
    let mut inv = prefix[n - 1].invert()?;

    // Compute inverses from right to left
    let mut result = vec![F::ZERO; n];
    for i in (1..n).rev() {
        result[i] = prefix[i - 1] * inv;
        inv = inv * elements[i];
    }
    result[0] = inv;

    Some(result)
}

/// Batch inversion in place, modifying the input slice.
/// Returns false if any element is zero.
pub fn batch_invert_inplace<F: Field>(elements: &mut [F]) -> bool {
    if let Some(inverted) = batch_invert(elements) {
        elements.copy_from_slice(&inverted);
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_id_roundtrip() {
        for b in [0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a] {
            let id = FieldId::from_byte(b).unwrap();
            assert_eq!(id.as_byte(), b);
        }
        assert!(FieldId::from_byte(0x00).is_none());
        assert!(FieldId::from_byte(0xff).is_none());
    }

    #[test]
    fn test_field_id_constants() {
        assert_eq!(Fp128::FIELD_ID, FieldId::Fp128);
        assert_eq!(GF2_128::FIELD_ID, FieldId::Gf2_128);
        assert_eq!(Fp256::FIELD_ID, FieldId::P256);
    }

    #[test]
    fn test_byte_lengths() {
        assert_eq!(Fp128::BYTE_LEN, 16);
        assert_eq!(GF2_128::BYTE_LEN, 16);
        assert_eq!(Fp256::BYTE_LEN, 32);
    }

    #[test]
    fn test_characteristics() {
        assert_eq!(GF2_128::CHARACTERISTIC, 2);
        // Fp128 and Fp256 have large primes that don't fit in u64
        assert_eq!(Fp128::CHARACTERISTIC, 0);
        assert_eq!(Fp256::CHARACTERISTIC, 0);
    }
}
