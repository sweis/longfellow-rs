//! Finite field arithmetic for the Longfellow ZK scheme.
//!
//! This module implements arithmetic over the field Fp128 = 2^128 - 2^108 + 1,
//! which is the primary field used in the Longfellow ZK protocol.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use zeroize::Zeroize;

/// The modulus p = 2^128 - 2^108 + 1
/// p = 340282042402384805036647824275747635201
/// In little-endian u64 limbs:
/// p & 0xFFFFFFFFFFFFFFFF = 1
/// p >> 64 = 0xfffff00000000000
const MODULUS_LO: u64 = 0x0000000000000001;
const MODULUS_HI: u64 = 0xfffff00000000000;

/// A field element in Fp128 = 2^128 - 2^108 + 1
///
/// Stored directly (not in Montgomery form) for simplicity and correctness.
#[derive(Clone, Copy, Debug, Zeroize)]
pub struct Fp128 {
    /// Low 64 bits of the value.
    lo: u64,
    /// High 64 bits of the value.
    hi: u64,
}

impl Default for Fp128 {
    fn default() -> Self {
        Self::ZERO
    }
}

impl ConstantTimeEq for Fp128 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.lo.ct_eq(&other.lo) & self.hi.ct_eq(&other.hi)
    }
}

impl ConditionallySelectable for Fp128 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            lo: u64::conditional_select(&a.lo, &b.lo, choice),
            hi: u64::conditional_select(&a.hi, &b.hi, choice),
        }
    }
}

impl PartialEq for Fp128 {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl Eq for Fp128 {}

impl Fp128 {
    /// The additive identity (zero).
    pub const ZERO: Self = Self { lo: 0, hi: 0 };

    /// The multiplicative identity (one).
    pub const ONE: Self = Self { lo: 1, hi: 0 };

    /// Two in the field (used for polynomial evaluation point P2).
    pub const TWO: Self = Self { lo: 2, hi: 0 };

    /// Create a new field element from a u64 value.
    pub fn from_u64(value: u64) -> Self {
        Self { lo: value, hi: 0 }
    }

    /// Create a field element from a u128 value.
    pub fn from_u128(value: u128) -> Self {
        let lo = value as u64;
        let hi = (value >> 64) as u64;
        Self { lo, hi }.reduce()
    }

    /// Create a field element from raw limbs.
    pub fn from_raw(limbs: [u64; 2]) -> Self {
        Self {
            lo: limbs[0],
            hi: limbs[1],
        }
        .reduce()
    }

    /// Create a field element from bytes (little-endian).
    pub fn from_bytes(bytes: &[u8; 16]) -> Self {
        let lo = u64::from_le_bytes(bytes[0..8].try_into().unwrap());
        let hi = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
        Self { lo, hi }.reduce()
    }

    /// Convert to bytes (little-endian).
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

    /// Reduce the value modulo p.
    fn reduce(self) -> Self {
        let mut result = self;

        // Check if we need to reduce (result >= p)
        while result.hi > MODULUS_HI
            || (result.hi == MODULUS_HI && result.lo >= MODULUS_LO)
        {
            // Subtract p
            let (lo, borrow) = result.lo.overflowing_sub(MODULUS_LO);
            let hi = result.hi.wrapping_sub(MODULUS_HI).wrapping_sub(borrow as u64);
            result = Self { lo, hi };
        }

        result
    }

    /// Compute the multiplicative inverse (a^{-1} mod p).
    /// Returns None if self is zero.
    pub fn invert(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // Using Fermat's little theorem: a^{-1} = a^{p-2} mod p
        // p - 2 = 2^128 - 2^108 - 1
        // p - 2 = 0xffffefffffffffffffffffffffffffff
        let p_minus_2 = Self {
            lo: 0xFFFFFFFFFFFFFFFF,
            hi: 0xffffefffffffffff,
        };

        Some(self.pow_vartime(&p_minus_2))
    }

    /// Square the field element.
    ///
    /// This is an optimized squaring operation that uses the fact that
    /// a^2 = a_hi^2 * 2^128 + 2*a_hi*a_lo * 2^64 + a_lo^2,
    /// requiring only 3 multiplications instead of 4 for general multiplication.
    pub fn square(&self) -> Self {
        // Compute 256-bit square as four 64-bit limbs
        // a^2 = (a_lo + a_hi * 2^64)^2
        //     = a_lo^2 + 2*a_lo*a_hi * 2^64 + a_hi^2 * 2^128
        let a0 = self.lo as u128;
        let a1 = self.hi as u128;

        // Partial products (3 instead of 4)
        let p00 = a0 * a0;          // a_lo^2, bits 0-127
        let p01 = a0 * a1;          // a_lo*a_hi, bits 64-191
        let p11 = a1 * a1;          // a_hi^2, bits 128-255

        // Combine into 256-bit result: [r0, r1, r2, r3]
        // result = r0 + r1*2^64 + r2*2^128 + r3*2^192
        // p01 appears twice (2*a_lo*a_hi), so we add it shifted left by 1 bit

        let r0 = p00 as u64;
        let carry0 = p00 >> 64;

        // p01 * 2 = 2*a_lo*a_hi at position 64
        // Check for overflow when doubling p01
        let p01_doubled = p01 << 1;
        let p01_doubled_overflow = (p01 >> 127) as u64; // 1 if p01 had bit 127 set

        let mid = carry0 + (p01_doubled as u64) as u128;
        let r1 = mid as u64;
        let carry1 = mid >> 64;

        let high = carry1 + (p01_doubled >> 64) + (p11 as u64) as u128 + (p01_doubled_overflow as u128) * (1u128 << 64);
        let r2 = high as u64;
        let carry2 = high >> 64;

        let r3 = (carry2 + (p11 >> 64)) as u64;

        // Reduce mod p = 2^128 - 2^108 + 1
        reduce_256_to_fp128(r0, r1, r2, r3)
    }

    /// Compute self^exp using square-and-multiply.
    fn pow_vartime(&self, exp: &Self) -> Self {
        let mut result = Self::ONE;
        let mut base = *self;

        // Process low 64 bits
        let mut e = exp.lo;
        for _ in 0..64 {
            if e & 1 == 1 {
                result = result * base;
            }
            base = base.square();
            e >>= 1;
        }

        // Process high 64 bits
        let mut e = exp.hi;
        for _ in 0..64 {
            if e & 1 == 1 {
                result = result * base;
            }
            base = base.square();
            e >>= 1;
        }

        result
    }

    /// Compute self^n for a u64 exponent.
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

    /// Check if the element is zero.
    pub fn is_zero(&self) -> bool {
        self.lo == 0 && self.hi == 0
    }

    /// Generate a random field element.
    pub fn random<R: rand::Rng>(rng: &mut R) -> Self {
        loop {
            let lo = rng.gen::<u64>();
            let hi = rng.gen::<u64>() & 0xffffff0000000000; // Mask to stay close to modulus

            let candidate = Self { lo, hi };
            if candidate.hi < MODULUS_HI
                || (candidate.hi == MODULUS_HI && candidate.lo < MODULUS_LO)
            {
                return candidate;
            }
        }
    }
}

impl Add for Fp128 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        // Add with carry tracking
        let (lo, carry1) = self.lo.overflowing_add(rhs.lo);
        let (hi, carry2) = self.hi.overflowing_add(rhs.hi);
        let (hi, carry3) = hi.overflowing_add(carry1 as u64);

        // If carry2 or carry3 is set, we have overflow past 2^128
        // We need to add (2^128 mod p) = 2^108 - 1 to the result
        let overflow = (carry2 as u64) | (carry3 as u64);

        if overflow == 0 {
            Self { lo, hi }.reduce()
        } else {
            // Add 2^108 - 1 for each overflow (can be at most 1 here since inputs < p)
            // 2^108 - 1 = 2^108 - 1 in limbs: lo = 0xFFFFFFFFFFFF (44 ones), hi = 0xFFF (12 ones at position 44)
            // Actually: 2^108 = 2^(64+44) = 2^44 in the hi part
            // 2^108 - 1 in 128 bits: low 108 bits are all 1
            // lo = 0xFFFFFFFFFFFFFFFF (all 64 bits), hi = 0x00000FFFFFFFFFFF (44 bits set)
            // Wait no: 2^108 - 1 has 108 ones
            // lo: 64 ones = 0xFFFFFFFFFFFFFFFF
            // hi: 44 ones = 0x00000FFFFFFFFFFF

            // 2^108 - 1 = {lo: 0xFFFFFFFFFFFFFFFF, hi: 0xFFFFFFFFFFF} but 0xFFFFFFFFFFF has 44 bits
            // 0xFFFFFFFFFFF = 2^44 - 1 = 17592186044415

            let correction_lo: u64 = 0xFFFFFFFFFFFFFFFF;
            let correction_hi: u64 = 0x00000FFFFFFFFFFF;

            let (lo2, c1) = lo.overflowing_add(correction_lo);
            let (hi2, c2) = hi.overflowing_add(correction_hi);
            let (hi2, c3) = hi2.overflowing_add(c1 as u64);

            // If there's still overflow, reduce again
            let overflow2 = (c2 as u64) | (c3 as u64);
            if overflow2 == 0 {
                Self { lo: lo2, hi: hi2 }.reduce()
            } else {
                // Apply correction again (very rare)
                let (lo3, c4) = lo2.overflowing_add(correction_lo);
                let (hi3, _) = hi2.overflowing_add(correction_hi);
                let (hi3, _) = hi3.overflowing_add(c4 as u64);
                Self { lo: lo3, hi: hi3 }.reduce()
            }
        }
    }
}

impl AddAssign for Fp128 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Fp128 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // If self < rhs, add p first
        if self.hi < rhs.hi || (self.hi == rhs.hi && self.lo < rhs.lo) {
            // self + p - rhs
            let (lo, carry) = self.lo.overflowing_add(MODULUS_LO);
            let hi = self.hi.wrapping_add(MODULUS_HI).wrapping_add(carry as u64);

            let (lo, borrow) = lo.overflowing_sub(rhs.lo);
            let hi = hi.wrapping_sub(rhs.hi).wrapping_sub(borrow as u64);

            Self { lo, hi }
        } else {
            // self - rhs
            let (lo, borrow) = self.lo.overflowing_sub(rhs.lo);
            let hi = self.hi.wrapping_sub(rhs.hi).wrapping_sub(borrow as u64);

            Self { lo, hi }
        }
    }
}

impl SubAssign for Fp128 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for Fp128 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.is_zero() {
            self
        } else {
            // p - self
            let (lo, borrow) = MODULUS_LO.overflowing_sub(self.lo);
            let hi = MODULUS_HI.wrapping_sub(self.hi).wrapping_sub(borrow as u64);
            Self { lo, hi }
        }
    }
}

impl Mul for Fp128 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // Compute 256-bit product as four 64-bit limbs
        // a * b = (a_lo + a_hi * 2^64) * (b_lo + b_hi * 2^64)
        let a0 = self.lo as u128;
        let a1 = self.hi as u128;
        let b0 = rhs.lo as u128;
        let b1 = rhs.hi as u128;

        // Partial products
        let p00 = a0 * b0;  // bits 0-127
        let p01 = a0 * b1;  // bits 64-191
        let p10 = a1 * b0;  // bits 64-191
        let p11 = a1 * b1;  // bits 128-255

        // Combine into 256-bit result: [r0, r1, r2, r3] where each is 64 bits
        // result = r0 + r1*2^64 + r2*2^128 + r3*2^192
        let r0 = p00 as u64;
        let carry0 = p00 >> 64;

        let mid = carry0 + (p01 as u64) as u128 + (p10 as u64) as u128;
        let r1 = mid as u64;
        let carry1 = mid >> 64;

        let high = carry1 + (p01 >> 64) + (p10 >> 64) + (p11 as u64) as u128;
        let r2 = high as u64;
        let carry2 = high >> 64;

        let r3 = (carry2 + (p11 >> 64)) as u64;

        // Reduce mod p = 2^128 - 2^108 + 1
        // Use the identity: 2^128 ≡ 2^108 - 1 (mod p)
        reduce_256_to_fp128(r0, r1, r2, r3)
    }
}

/// Reduce a 256-bit value (represented as 4 x 64-bit limbs) modulo p.
/// value = r0 + r1*2^64 + r2*2^128 + r3*2^192
fn reduce_256_to_fp128(r0: u64, r1: u64, r2: u64, r3: u64) -> Fp128 {
    // p = 2^128 - 2^108 + 1
    // 2^128 ≡ 2^108 - 1 (mod p)
    //
    // Strategy: reduce r3*2^192 and r2*2^128 using the identity,
    // accumulating into a 192-bit intermediate, then reduce again.

    // Start with the low 128 bits
    let mut lo: u128 = (r0 as u128) | ((r1 as u128) << 64);
    let mut hi: u128 = 0;  // bits 128+

    // Add r2 * 2^128 ≡ r2 * (2^108 - 1) (mod p)
    // = r2 * 2^108 - r2
    if r2 > 0 {
        let r2_128 = r2 as u128;

        // r2 * 2^108: this is at most 64+108 = 172 bits
        // Split: (r2 * 2^108) = (r2 << 108)
        // low 128 bits and high bits
        let shifted_lo = r2_128 << 108;  // only low 20 bits of r2 contribute to low 128
        let shifted_hi = r2_128 >> 20;   // bits that overflow into position 128+

        let (new_lo, c1) = lo.overflowing_add(shifted_lo);
        lo = new_lo;
        hi = hi.wrapping_add(shifted_hi).wrapping_add(c1 as u128);

        // Subtract r2
        let (new_lo, b1) = lo.overflowing_sub(r2_128);
        lo = new_lo;
        if b1 {
            hi = hi.wrapping_sub(1);
        }
    }

    // Add r3 * 2^192
    // 2^192 = 2^64 * 2^128 ≡ 2^64 * (2^108 - 1) = 2^172 - 2^64 (mod p)
    // 2^172 = 2^44 * 2^128 ≡ 2^44 * (2^108 - 1) = 2^152 - 2^44 (mod p)
    // 2^152 = 2^24 * 2^128 ≡ 2^24 * (2^108 - 1) = 2^132 - 2^24 (mod p)
    // 2^132 = 2^4 * 2^128 ≡ 2^4 * (2^108 - 1) = 2^112 - 2^4 (mod p)
    // So: 2^152 ≡ 2^112 - 2^4 - 2^24 (mod p)
    // So: 2^172 ≡ 2^112 - 2^4 - 2^24 - 2^44 (mod p)
    // So: 2^192 ≡ 2^112 - 2^4 - 2^24 - 2^44 - 2^64 (mod p)
    if r3 > 0 {
        let r3_128 = r3 as u128;

        // Add r3 * 2^112
        let add_term = r3_128 << 112;  // at most 64+112 = 176 bits
        let add_lo = add_term;         // low 128 bits (wraps)
        let add_hi = r3_128 >> 16;     // bits that go above 128

        let (new_lo, c) = lo.overflowing_add(add_lo);
        lo = new_lo;
        hi = hi.wrapping_add(add_hi).wrapping_add(c as u128);

        // Subtract r3 * (2^64 + 2^44 + 2^24 + 2^4)
        // These all fit in 128 bits for reasonable r3
        let sub_term = r3_128
            .wrapping_mul(
                (1u128 << 64)
                    .wrapping_add(1u128 << 44)
                    .wrapping_add(1u128 << 24)
                    .wrapping_add(1u128 << 4),
            );

        let (new_lo, b) = lo.overflowing_sub(sub_term);
        lo = new_lo;
        if b {
            hi = hi.wrapping_sub(1);
        }
    }

    // Now we have a value in (lo, hi) where hi represents bits 128+
    // Keep reducing hi using: hi * 2^128 ≡ hi * (2^108 - 1) (mod p)
    let mut iterations = 0;
    while hi != 0 && iterations < 10 {
        iterations += 1;
        let h = hi;
        hi = 0;

        // Add h * 2^108
        let h_shifted_lo = h << 108;
        let h_shifted_hi = h >> 20;

        let (new_lo, c) = lo.overflowing_add(h_shifted_lo);
        lo = new_lo;
        hi = hi.wrapping_add(h_shifted_hi).wrapping_add(c as u128);

        // Subtract h
        let (new_lo, b) = lo.overflowing_sub(h);
        lo = new_lo;
        if b {
            hi = hi.wrapping_sub(1);
        }
    }

    // If hi wrapped around to a large value (negative), add p back
    if hi > 0x1000000 {
        // This means we went negative - add p
        let p = ((MODULUS_HI as u128) << 64) | (MODULUS_LO as u128);
        lo = lo.wrapping_add(p);
    }

    // Final reduction: ensure lo < p
    let p = ((MODULUS_HI as u128) << 64) | (MODULUS_LO as u128);
    while lo >= p {
        lo = lo.wrapping_sub(p);
    }

    Fp128 {
        lo: lo as u64,
        hi: (lo >> 64) as u64,
    }
}

impl MulAssign for Fp128 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl std::iter::Sum for Fp128 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |acc, x| acc + x)
    }
}

impl std::iter::Product for Fp128 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ONE, |acc, x| acc * x)
    }
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
        // result[i] = (a_0 * ... * a_{i-1})^{-1} * (a_{i+1} * ... * a_{n-1})^{-1}
        //           = prefix[i-1] * suffix_inv
        result[i] = prefix[i - 1] * inv;
        // Update suffix inverse: suffix_inv = a_i * suffix_inv_{i+1}
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

    /// The evaluation point P2 for polynomial representation.
    /// For prime fields with characteristic > 2, this is 2.
    const P2: Self;

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

    /// Convert to bytes.
    fn to_bytes(&self) -> Vec<u8>;

    /// Create from bytes.
    fn from_bytes(bytes: &[u8]) -> Self;
}

impl Field for Fp128 {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;
    const P2: Self = Self::TWO;

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
        assert!(Fp128::ZERO.is_zero());
        assert!(!Fp128::ONE.is_zero());
        assert_eq!(Fp128::ZERO + Fp128::ONE, Fp128::ONE);
        assert_eq!(Fp128::ONE - Fp128::ONE, Fp128::ZERO);
    }

    #[test]
    fn test_addition() {
        let a = Fp128::from_u64(5);
        let b = Fp128::from_u64(7);
        let c = Fp128::from_u64(12);
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_subtraction() {
        let a = Fp128::from_u64(10);
        let b = Fp128::from_u64(7);
        let c = Fp128::from_u64(3);
        assert_eq!(a - b, c);
    }

    #[test]
    fn test_multiplication() {
        let a = Fp128::from_u64(5);
        let b = Fp128::from_u64(7);
        let c = Fp128::from_u64(35);
        assert_eq!(a * b, c);
    }

    #[test]
    fn test_negation() {
        let a = Fp128::from_u64(5);
        assert_eq!(a + (-a), Fp128::ZERO);
    }

    #[test]
    fn test_inversion() {
        let a = Fp128::from_u64(7);
        let a_inv = a.invert().unwrap();
        assert_eq!(a * a_inv, Fp128::ONE);
    }

    #[test]
    fn test_two() {
        let two = Fp128::from_u64(2);
        assert_eq!(two, Fp128::ONE + Fp128::ONE);
    }

    #[test]
    fn test_batch_invert() {
        // Test batch inversion
        let elements = vec![
            Fp128::from_u64(2),
            Fp128::from_u64(3),
            Fp128::from_u64(5),
            Fp128::from_u64(7),
            Fp128::from_u64(11),
        ];

        let inverted = super::batch_invert(&elements).unwrap();

        // Verify each inverse
        for (a, a_inv) in elements.iter().zip(inverted.iter()) {
            assert_eq!(*a * *a_inv, Fp128::ONE);
        }
    }

    #[test]
    fn test_batch_invert_single() {
        let elements = vec![Fp128::from_u64(7)];
        let inverted = super::batch_invert(&elements).unwrap();
        assert_eq!(elements[0] * inverted[0], Fp128::ONE);
    }

    #[test]
    fn test_batch_invert_empty() {
        let elements: Vec<Fp128> = vec![];
        let inverted = super::batch_invert(&elements).unwrap();
        assert!(inverted.is_empty());
    }

    #[test]
    fn test_batch_invert_zero() {
        let elements = vec![Fp128::from_u64(2), Fp128::ZERO, Fp128::from_u64(3)];
        assert!(super::batch_invert(&elements).is_none());
    }

    #[test]
    fn test_optimized_squaring() {
        // Test that optimized squaring matches multiplication
        let test_values = [
            Fp128::from_u64(0),
            Fp128::from_u64(1),
            Fp128::from_u64(2),
            Fp128::from_u64(7),
            Fp128::from_u64(12345),
            Fp128::from_u64(0xFFFFFFFFFFFFFFFF),
            Fp128::from_raw([0xFFFFFFFFFFFFFFFF, 0x00000F0000000000]),
            Fp128::from_raw([0x123456789ABCDEF0, 0x00000FEDCBA98765]),
        ];

        for &a in &test_values {
            let squared = a.square();
            let multiplied = a * a;
            assert_eq!(squared, multiplied, "Squaring mismatch for {:?}", a);
        }
    }

    #[test]
    fn test_squaring_random() {
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        for _ in 0..100 {
            let a = Fp128::random(&mut rng);
            let squared = a.square();
            let multiplied = a * a;
            assert_eq!(squared, multiplied);
        }
    }
}
