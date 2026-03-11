//! Binary extension field GF(2^128).
//!
//! This is FieldID 0x04 in the IETF draft-google-cfrg-libzk specification.
//!
//! The field is defined as GF(2)[x] / Q(x) where:
//!   Q(x) = x^128 + x^7 + x^2 + x + 1
//!
//! With this choice, x is a generator of the multiplicative group.
//!
//! # Key properties for ZK:
//! - Characteristic 2: a + a = 0 for all a (addition is XOR)
//! - Negation is the identity: -a = a
//! - P2 = X (the polynomial variable, bit representation 0b10)
//!
//! # Multiplication
//! Uses carryless (polynomial) multiplication followed by reduction
//! modulo Q(x). The reduction uses the identity x^128 ≡ x^7 + x^2 + x + 1.

use super::{Field, FieldId};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use zeroize::Zeroize;

/// The reduction polynomial tail: x^7 + x^2 + x + 1 = 0b10000111 = 0x87.
/// This represents Q(x) - x^128, i.e., the low-degree part of the
/// irreducible polynomial.
const REDUCTION_POLY: u64 = 0x87;

/// A field element in GF(2^128) = GF(2)[x] / (x^128 + x^7 + x^2 + x + 1).
///
/// Elements are stored as 128-bit polynomials over GF(2), packed into
/// two u64 limbs in little-endian bit order (lo = bits 0..63, hi = bits 64..127).
#[derive(Clone, Copy, Debug, Zeroize)]
#[allow(non_camel_case_types)]
pub struct GF2_128 {
    lo: u64,
    hi: u64,
}

impl Default for GF2_128 {
    fn default() -> Self {
        Self::ZERO
    }
}

impl ConstantTimeEq for GF2_128 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.lo.ct_eq(&other.lo) & self.hi.ct_eq(&other.hi)
    }
}

impl ConditionallySelectable for GF2_128 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            lo: u64::conditional_select(&a.lo, &b.lo, choice),
            hi: u64::conditional_select(&a.hi, &b.hi, choice),
        }
    }
}

impl PartialEq for GF2_128 {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl Eq for GF2_128 {}

impl GF2_128 {
    /// The additive identity (zero polynomial).
    pub const ZERO: Self = Self { lo: 0, hi: 0 };

    /// The multiplicative identity (constant polynomial 1).
    pub const ONE: Self = Self { lo: 1, hi: 0 };

    /// The element X (the polynomial variable x).
    /// This is the generator and also the P2 evaluation point for sumcheck.
    pub const X: Self = Self { lo: 2, hi: 0 };

    /// x^{-1} = x^127 + x^6 + x + 1.
    /// Derived from: x * x^{-1} = 1, and x^128 = x^7 + x^2 + x + 1.
    /// So x^{-1} = (x^7 + x^2 + x + 1) / x... but since we're in GF(2),
    /// x * (x^127 + x^6 + x + 1) = x^128 + x^7 + x^2 + x
    ///   = (x^7 + x^2 + x + 1) + x^7 + x^2 + x = 1. ✓
    pub const X_INV: Self = Self {
        lo: 0x43, // 0b1000011 = x^6 + x + 1
        hi: 0x8000000000000000, // x^127
    };

    /// Create a field element from a u64 value (treated as a degree < 64 polynomial).
    pub fn from_u64(value: u64) -> Self {
        Self { lo: value, hi: 0 }
    }

    /// Create a field element from a u128 value.
    pub fn from_u128(value: u128) -> Self {
        Self {
            lo: value as u64,
            hi: (value >> 64) as u64,
        }
    }

    /// Create from raw limbs (no reduction needed since any 128-bit value is valid).
    pub fn from_raw(limbs: [u64; 2]) -> Self {
        Self {
            lo: limbs[0],
            hi: limbs[1],
        }
    }

    /// Create from 16 little-endian bytes.
    pub fn from_bytes(bytes: &[u8; 16]) -> Self {
        let lo = u64::from_le_bytes(bytes[0..8].try_into().unwrap());
        let hi = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
        Self { lo, hi }
    }

    /// Convert to 16 little-endian bytes.
    pub fn to_bytes(&self) -> [u8; 16] {
        let mut bytes = [0u8; 16];
        bytes[0..8].copy_from_slice(&self.lo.to_le_bytes());
        bytes[8..16].copy_from_slice(&self.hi.to_le_bytes());
        bytes
    }

    /// Get the value as a u128.
    pub fn to_u128(&self) -> u128 {
        (self.hi as u128) << 64 | (self.lo as u128)
    }

    /// Check if the element is zero.
    pub fn is_zero(&self) -> bool {
        self.lo == 0 && self.hi == 0
    }

    /// Square the element.
    ///
    /// In GF(2^k), squaring is a linear operation (Frobenius endomorphism).
    /// (a + b)^2 = a^2 + b^2 (no cross term since 2ab = 0 in char 2).
    /// This means we just need to "spread out" the bits: each bit i
    /// moves to position 2i, then reduce.
    pub fn square(&self) -> Self {
        // Spread bits: bit i -> position 2i
        // This produces a 256-bit result that needs reduction.
        let (r0, r1) = spread_bits_64(self.lo); // bits 0..63 -> 0..127
        let (r2, r3) = spread_bits_64(self.hi); // bits 64..127 -> 128..255

        reduce_256(r0, r1, r2, r3)
    }

    /// Compute the multiplicative inverse using Fermat's little theorem.
    ///
    /// In GF(2^128), the multiplicative group has order 2^128 - 1, so
    /// a^{-1} = a^{2^128 - 2}.
    ///
    /// We use an addition chain optimized for the exponent 2^128 - 2
    /// which has the binary representation of 127 ones followed by a zero:
    /// 0b111...1110 (127 ones, then 0).
    pub fn invert(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // Use Itoh-Tsujii inversion: a^(-1) = a^(2^128 - 2)
        // 2^128 - 2 = 2 * (2^127 - 1)
        //
        // First compute a^(2^127 - 1) using an addition chain,
        // then square once.
        //
        // 2^127 - 1 = (2^64 - 1) * 2^63 + (2^63 - 1)
        //           ... use standard doubling chain.
        //
        // We compute a^(2^k - 1) iteratively using:
        // a^(2^(2k) - 1) = (a^(2^k - 1))^(2^k) * a^(2^k - 1)

        let a = *self;

        // a^(2^1 - 1) = a
        let a1 = a;

        // a^(2^2 - 1) = a^3 = a^2 * a
        let a2 = a1.square() * a1;

        // a^(2^4 - 1) = (a^(2^2-1))^(2^2) * a^(2^2-1)
        let mut t = a2;
        for _ in 0..2 { t = t.square(); }
        let a4 = t * a2;

        // a^(2^8 - 1)
        let mut t = a4;
        for _ in 0..4 { t = t.square(); }
        let a8 = t * a4;

        // a^(2^16 - 1)
        let mut t = a8;
        for _ in 0..8 { t = t.square(); }
        let a16 = t * a8;

        // a^(2^32 - 1)
        let mut t = a16;
        for _ in 0..16 { t = t.square(); }
        let a32 = t * a16;

        // a^(2^64 - 1)
        let mut t = a32;
        for _ in 0..32 { t = t.square(); }
        let a64 = t * a32;

        // a^(2^127 - 1) = (a^(2^64 - 1))^(2^63) * a^(2^63 - 1)
        // We need a^(2^63 - 1) first.
        // a^(2^63 - 1) = (a^(2^32 - 1))^(2^31) * a^(2^31 - 1)
        // a^(2^31 - 1) = (a^(2^16 - 1))^(2^15) * a^(2^15 - 1)
        // a^(2^15 - 1) = (a^(2^8 - 1))^(2^7) * a^(2^7 - 1)
        // a^(2^7 - 1) = (a^(2^4 - 1))^(2^3) * a^(2^3 - 1)
        // a^(2^3 - 1) = (a^(2^2 - 1))^(2^1) * a^(2^1 - 1)

        let mut t = a2;
        t = t.square();
        let a3 = t * a1; // a^(2^3 - 1)

        let mut t = a4;
        for _ in 0..3 { t = t.square(); }
        let a7 = t * a3; // a^(2^7 - 1)

        let mut t = a8;
        for _ in 0..7 { t = t.square(); }
        let a15 = t * a7; // a^(2^15 - 1)

        let mut t = a16;
        for _ in 0..15 { t = t.square(); }
        let a31 = t * a15; // a^(2^31 - 1)

        let mut t = a32;
        for _ in 0..31 { t = t.square(); }
        let a63 = t * a31; // a^(2^63 - 1)

        let mut t = a64;
        for _ in 0..63 { t = t.square(); }
        let a127 = t * a63; // a^(2^127 - 1)

        // Finally: a^(-1) = a^(2^128 - 2) = (a^(2^127 - 1))^2
        Some(a127.square())
    }

    /// Compute self^n using square-and-multiply.
    pub fn pow(&self, n: u64) -> Self {
        let mut result = Self::ONE;
        let mut base = *self;
        let mut exp = n;

        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base.square();
            exp >>= 1;
        }

        result
    }

    /// Generate a random field element.
    /// All 128-bit values are valid field elements, so this is simple.
    pub fn random<R: rand::Rng>(rng: &mut R) -> Self {
        Self {
            lo: rng.gen::<u64>(),
            hi: rng.gen::<u64>(),
        }
    }
}

/// Spread 64 bits so that bit i goes to position 2i.
/// Returns (low 64 bits, high 64 bits) of the 128-bit result.
///
/// Used for squaring in GF(2^k): each coefficient x^i in the input
/// becomes x^(2i) in the square.
fn spread_bits_64(x: u64) -> (u64, u64) {
    // Spread the low 32 bits into the low 64 bits of output
    // Spread the high 32 bits into the high 64 bits of output
    (spread_bits_32(x as u32), spread_bits_32((x >> 32) as u32))
}

/// Spread 32 bits so that bit i goes to position 2i in a 64-bit result.
fn spread_bits_32(x: u32) -> u64 {
    let mut x = x as u64;
    // Classic bit-interleaving trick (Morton encoding)
    x = (x | (x << 16)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x << 8))  & 0x00FF_00FF_00FF_00FF;
    x = (x | (x << 4))  & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x << 2))  & 0x3333_3333_3333_3333;
    x = (x | (x << 1))  & 0x5555_5555_5555_5555;
    x
}

/// Portable 64x64 -> 128 carryless (polynomial) multiplication.
///
/// Computes the product of two polynomials over GF(2) with degrees < 64.
/// The result has degree < 127, fitting in 128 bits.
fn clmul64(a: u64, b: u64) -> (u64, u64) {
    // Simple bit-by-bit implementation.
    // For performance-critical code, this should use CLMUL intrinsics
    // (PCLMULQDQ on x86, PMULL on ARM).
    let mut lo: u64 = 0;
    let mut hi: u64 = 0;

    for i in 0..64 {
        if (b >> i) & 1 == 1 {
            // Add a << i (XOR in GF(2))
            if i == 0 {
                lo ^= a;
            } else {
                lo ^= a << i;
                hi ^= a >> (64 - i);
            }
        }
    }

    (lo, hi)
}

/// Reduce a 256-bit polynomial modulo Q(x) = x^128 + x^7 + x^2 + x + 1.
///
/// Input: 256-bit value as four 64-bit limbs [r0, r1, r2, r3]
/// representing r0 + r1*x^64 + r2*x^128 + r3*x^192.
///
/// Uses the identity: x^128 ≡ x^7 + x^2 + x + 1 (mod Q(x)).
fn reduce_256(r0: u64, r1: u64, r2: u64, r3: u64) -> GF2_128 {
    // We need to fold r2*x^128 + r3*x^192 back into the low 128 bits.
    // x^128 ≡ x^7 + x^2 + x + 1 = REDUCTION_POLY (0x87)
    // x^192 = x^64 * x^128 ≡ x^64 * (x^7 + x^2 + x + 1)
    //
    // Strategy: for each high limb, multiply by the reduction polynomial
    // and XOR into the low limbs.

    // First, fold r3*x^192:
    // r3*x^192 ≡ r3 * x^64 * (x^7 + x^2 + x + 1)
    //          = (r3 * 0x87) shifted to position x^64
    // The product r3 * 0x87 is at most 64 + 8 = 72 bits.
    let (p3_lo, p3_hi) = clmul64(r3, REDUCTION_POLY);
    // p3_lo goes to positions 64..127 (limb 1)
    // p3_hi goes to positions 128..191 (limb 2) - needs further reduction!

    let r1 = r1 ^ p3_lo;
    let r2 = r2 ^ p3_hi;

    // Now fold r2*x^128:
    // r2*x^128 ≡ r2 * (x^7 + x^2 + x + 1) = r2 * 0x87
    // The product is at most 64 + 8 = 72 bits.
    let (p2_lo, p2_hi) = clmul64(r2, REDUCTION_POLY);
    // p2_lo goes to positions 0..63 (limb 0)
    // p2_hi goes to positions 64..71 (limb 1, low 8 bits) - fits entirely!

    let r0 = r0 ^ p2_lo;
    let r1 = r1 ^ p2_hi;

    // But wait: p3_hi could have contributed up to 8 bits to r2,
    // which when multiplied by 0x87 gives at most 8+8 = 16 bits,
    // all of which fit in the low limbs. So no further reduction needed.

    GF2_128 { lo: r0, hi: r1 }
}

// === Arithmetic operations ===

impl Add for GF2_128 {
    type Output = Self;

    /// Addition in GF(2^k) is XOR (componentwise addition mod 2).
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            lo: self.lo ^ rhs.lo,
            hi: self.hi ^ rhs.hi,
        }
    }
}

impl AddAssign for GF2_128 {
    fn add_assign(&mut self, rhs: Self) {
        self.lo ^= rhs.lo;
        self.hi ^= rhs.hi;
    }
}

impl Sub for GF2_128 {
    type Output = Self;

    /// Subtraction in GF(2^k) is the same as addition (since -1 = 1).
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Self) -> Self::Output {
        self + rhs
    }
}

impl SubAssign for GF2_128 {
    #[allow(clippy::suspicious_op_assign_impl)]
    fn sub_assign(&mut self, rhs: Self) {
        *self += rhs;
    }
}

impl Neg for GF2_128 {
    type Output = Self;

    /// Negation in GF(2^k) is the identity (since -1 = 1).
    fn neg(self) -> Self::Output {
        self
    }
}

impl Mul for GF2_128 {
    type Output = Self;

    /// Multiplication: carryless multiply then reduce modulo Q(x).
    fn mul(self, rhs: Self) -> Self::Output {
        // Compute the full 256-bit carryless product
        // (a_lo + a_hi*x^64) * (b_lo + b_hi*x^64)
        // = a_lo*b_lo + (a_lo*b_hi + a_hi*b_lo)*x^64 + a_hi*b_hi*x^128

        let (p00_lo, p00_hi) = clmul64(self.lo, rhs.lo);
        let (p01_lo, p01_hi) = clmul64(self.lo, rhs.hi);
        let (p10_lo, p10_hi) = clmul64(self.hi, rhs.lo);
        let (p11_lo, p11_hi) = clmul64(self.hi, rhs.hi);

        // Combine (all additions are XOR):
        // r0 = p00_lo
        // r1 = p00_hi ^ p01_lo ^ p10_lo
        // r2 = p01_hi ^ p10_hi ^ p11_lo
        // r3 = p11_hi
        let r0 = p00_lo;
        let r1 = p00_hi ^ p01_lo ^ p10_lo;
        let r2 = p01_hi ^ p10_hi ^ p11_lo;
        let r3 = p11_hi;

        reduce_256(r0, r1, r2, r3)
    }
}

impl MulAssign for GF2_128 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl std::iter::Sum for GF2_128 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |acc, x| acc + x)
    }
}

impl std::iter::Product for GF2_128 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ONE, |acc, x| acc * x)
    }
}

impl Field for GF2_128 {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;
    /// P2 = X per the spec: for GF(2^128), the third sumcheck evaluation
    /// point is the polynomial variable X.
    const P2: Self = Self::X;
    const FIELD_ID: FieldId = FieldId::Gf2_128;
    const BYTE_LEN: usize = 16;
    const CHARACTERISTIC: u64 = 2;

    fn from_u64(value: u64) -> Self {
        Self::from_u64(value)
    }

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn invert(&self) -> Option<Self> {
        self.invert()
    }

    fn square(&self) -> Self {
        self.square()
    }

    fn random<R: rand::Rng>(rng: &mut R) -> Self {
        Self::random(rng)
    }

    fn to_bytes(&self) -> Vec<u8> {
        self.to_bytes().to_vec()
    }

    fn from_bytes(bytes: &[u8]) -> Self {
        let mut arr = [0u8; 16];
        let len = bytes.len().min(16);
        arr[..len].copy_from_slice(&bytes[..len]);
        Self::from_bytes(&arr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_one() {
        assert!(GF2_128::ZERO.is_zero());
        assert!(!GF2_128::ONE.is_zero());
        assert_eq!(GF2_128::ZERO + GF2_128::ONE, GF2_128::ONE);
        assert_eq!(GF2_128::ONE + GF2_128::ONE, GF2_128::ZERO); // char 2!
    }

    #[test]
    fn test_addition_is_xor() {
        let a = GF2_128::from_u64(0b1010);
        let b = GF2_128::from_u64(0b0110);
        assert_eq!(a + b, GF2_128::from_u64(0b1100));
    }

    #[test]
    fn test_subtraction_equals_addition() {
        let a = GF2_128::from_u64(123);
        let b = GF2_128::from_u64(456);
        assert_eq!(a - b, a + b);
    }

    #[test]
    fn test_negation_is_identity() {
        let a = GF2_128::from_u64(12345);
        assert_eq!(-a, a);
        assert_eq!(a + (-a), GF2_128::ZERO);
    }

    #[test]
    fn test_multiplication_small() {
        // x * x = x^2 (no reduction needed for small values)
        let x = GF2_128::X; // = 0b10
        let x2 = x * x;
        assert_eq!(x2, GF2_128::from_u64(0b100)); // x^2

        // x^2 * x = x^3
        let x3 = x2 * x;
        assert_eq!(x3, GF2_128::from_u64(0b1000));
    }

    #[test]
    fn test_multiplication_by_one() {
        let a = GF2_128::from_raw([0x123456789abcdef0, 0xfedcba9876543210]);
        assert_eq!(a * GF2_128::ONE, a);
        assert_eq!(GF2_128::ONE * a, a);
    }

    #[test]
    fn test_multiplication_by_zero() {
        let a = GF2_128::from_raw([0x123456789abcdef0, 0xfedcba9876543210]);
        assert_eq!(a * GF2_128::ZERO, GF2_128::ZERO);
    }

    #[test]
    fn test_multiplication_commutative() {
        let a = GF2_128::from_raw([0x123456789abcdef0, 0xfedcba9876543210]);
        let b = GF2_128::from_raw([0xdeadbeefcafebabe, 0x0011223344556677]);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_multiplication_associative() {
        let a = GF2_128::from_u64(7);
        let b = GF2_128::from_u64(13);
        let c = GF2_128::from_u64(17);
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn test_distributive() {
        let a = GF2_128::from_u64(5);
        let b = GF2_128::from_u64(7);
        let c = GF2_128::from_u64(11);
        assert_eq!(a * (b + c), a * b + a * c);
    }

    #[test]
    fn test_reduction() {
        // x^128 should reduce to x^7 + x^2 + x + 1 = 0x87
        // We can test this by multiplying x^64 * x^64
        let x64 = GF2_128::from_raw([0, 1]); // x^64
        let x128_reduced = x64 * x64;
        assert_eq!(x128_reduced, GF2_128::from_u64(0x87));
    }

    #[test]
    fn test_x_inverse() {
        // X * X_INV should equal 1
        assert_eq!(GF2_128::X * GF2_128::X_INV, GF2_128::ONE);
    }

    #[test]
    fn test_inversion() {
        // Test several values
        for v in [1u64, 2, 3, 5, 7, 0x87, 0xdeadbeef, u64::MAX] {
            let a = GF2_128::from_u64(v);
            let a_inv = a.invert().expect("non-zero element should be invertible");
            assert_eq!(a * a_inv, GF2_128::ONE, "Failed for v={}", v);
        }
    }

    #[test]
    fn test_inversion_large() {
        let a = GF2_128::from_raw([0x123456789abcdef0, 0xfedcba9876543210]);
        let a_inv = a.invert().unwrap();
        assert_eq!(a * a_inv, GF2_128::ONE);
        // Double inversion
        assert_eq!(a_inv.invert().unwrap(), a);
    }

    #[test]
    fn test_inversion_zero() {
        assert!(GF2_128::ZERO.invert().is_none());
    }

    #[test]
    fn test_squaring_matches_multiplication() {
        let test_values = [
            GF2_128::ZERO,
            GF2_128::ONE,
            GF2_128::X,
            GF2_128::from_u64(7),
            GF2_128::from_u64(0x87),
            GF2_128::from_raw([0xFFFFFFFFFFFFFFFF, 0]),
            GF2_128::from_raw([0, 0xFFFFFFFFFFFFFFFF]),
            GF2_128::from_raw([0x123456789ABCDEF0, 0xFEDCBA9876543210]),
        ];

        for &a in &test_values {
            assert_eq!(a.square(), a * a, "Squaring mismatch for {:?}", a);
        }
    }

    #[test]
    fn test_frobenius() {
        // In GF(2^k), (a + b)^2 = a^2 + b^2 (Frobenius is additive)
        let a = GF2_128::from_u64(0x5555);
        let b = GF2_128::from_u64(0xaaaa);
        assert_eq!((a + b).square(), a.square() + b.square());
    }

    #[test]
    fn test_spread_bits() {
        // 0b1011 should spread to 0b01000101
        assert_eq!(spread_bits_32(0b1011), 0b01000101);
        // All ones should spread to 0x5555...
        assert_eq!(spread_bits_32(0xFFFFFFFF), 0x5555_5555_5555_5555);
    }

    #[test]
    fn test_clmul_basic() {
        // 0b11 * 0b11 = x^2 + 1 (not x^2 + 2x + 1, no carries!)
        // = (x+1)*(x+1) = x^2 + 2x + 1 = x^2 + 1 in GF(2)
        let (lo, hi) = clmul64(0b11, 0b11);
        assert_eq!(lo, 0b101);
        assert_eq!(hi, 0);
    }

    #[test]
    fn test_serialization() {
        let a = GF2_128::from_raw([0x123456789abcdef0, 0xfedcba9876543210]);
        let bytes = a.to_bytes();
        assert_eq!(bytes.len(), 16);
        let recovered = GF2_128::from_bytes(&bytes);
        assert_eq!(a, recovered);
    }

    #[test]
    fn test_random() {
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let a = GF2_128::random(&mut rng);
        let b = GF2_128::random(&mut rng);
        assert_ne!(a, b); // Extremely unlikely to be equal
    }

    #[test]
    fn test_p2_constant() {
        // P2 should be X, not the same as ONE+ONE (which is ZERO in char 2!)
        assert_eq!(<GF2_128 as Field>::P2, GF2_128::X);
        assert_ne!(<GF2_128 as Field>::P2, GF2_128::ZERO);
        assert_ne!(<GF2_128 as Field>::P2, GF2_128::ONE);
        // Sanity: in char 2, 1+1 = 0, so P2 can't be "two"
        assert_eq!(GF2_128::ONE + GF2_128::ONE, GF2_128::ZERO);
    }

    #[test]
    fn test_pow() {
        let x = GF2_128::X;
        // x^7 * x^121 = x^128 = 0x87
        assert_eq!(x.pow(128), GF2_128::from_u64(0x87));
        // x^0 = 1
        assert_eq!(x.pow(0), GF2_128::ONE);
        // x^1 = x
        assert_eq!(x.pow(1), x);
    }

    #[test]
    fn test_batch_invert() {
        let elements = vec![
            GF2_128::from_u64(2),
            GF2_128::from_u64(3),
            GF2_128::from_u64(5),
            GF2_128::from_raw([0xdeadbeef, 0xcafebabe]),
        ];

        let inverted = crate::field::batch_invert(&elements).unwrap();
        for (a, a_inv) in elements.iter().zip(inverted.iter()) {
            assert_eq!(*a * *a_inv, GF2_128::ONE);
        }
    }
}
