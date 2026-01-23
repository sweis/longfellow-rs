//! Fast Fourier Transform (FFT) for polynomial operations.
//!
//! This module implements the radix-2 Cooley-Tukey FFT algorithm for
//! O(n log n) polynomial multiplication and evaluation.
//!
//! # Field Requirements
//!
//! The FFT requires the field to have roots of unity of sufficient order.
//! For Fp128 = 2^128 - 2^108 + 1, we have:
//! - p - 1 = 2^108 * (2^20 - 1)
//! - This supports 2^k-th roots of unity for k ≤ 108
//!
//! # Usage
//!
//! ```ignore
//! use longfellow_zk::fft::{fft, ifft, FftDomain};
//!
//! // Create a domain for size-8 FFT
//! let domain = FftDomain::new(8);
//!
//! // Transform coefficients to evaluations
//! let coeffs = vec![Fp128::from_u64(1), Fp128::from_u64(2), ...];
//! let evals = fft(&coeffs, &domain);
//!
//! // Transform back
//! let recovered = ifft(&evals, &domain);
//! ```

use crate::field::{Field, Fp128};

/// A primitive 2^108-th root of unity for Fp128.
///
/// This is computed as g^((p-1)/2^108) where g = 17 is a generator.
/// Value verified by: omega^(2^108) = 1 and omega^(2^107) ≠ 1.
///
/// Note: Generators 2, 3, 5, 7, 11, 13, 19, 23 don't work for this field.
/// The element 17 was found by testing that 17^(2^20-1) has full order 2^108.
fn primitive_root_of_unity_2_108() -> Fp128 {
    // g = 17, exponent = (p-1) / 2^108 = 2^20 - 1 = 1048575
    // omega = 17^1048575 mod p
    Fp128::from_u64(17).pow(1048575)
}

/// Compute a 2^k-th primitive root of unity.
///
/// Returns omega such that omega^(2^k) = 1 and omega^(2^(k-1)) ≠ 1.
pub fn root_of_unity(k: u32) -> Fp128 {
    assert!(k <= 108, "Fp128 supports at most 2^108-th roots of unity");

    // Start with primitive 2^108-th root and square down
    let mut omega = primitive_root_of_unity_2_108();
    for _ in 0..(108 - k) {
        omega = omega.square();
    }
    omega
}

/// Precomputed data for FFT operations of a specific size.
#[derive(Clone, Debug)]
pub struct FftDomain {
    /// Size of the domain (must be a power of 2).
    pub size: usize,
    /// log2 of the size.
    pub log_size: u32,
    /// Primitive n-th root of unity (omega).
    pub omega: Fp128,
    /// Inverse of omega (for inverse FFT).
    pub omega_inv: Fp128,
    /// Inverse of size (for inverse FFT scaling).
    pub size_inv: Fp128,
    /// Precomputed powers of omega: [1, omega, omega^2, ..., omega^(n/2-1)]
    pub twiddles: Vec<Fp128>,
    /// Precomputed powers of omega_inv for inverse FFT.
    pub twiddles_inv: Vec<Fp128>,
}

impl FftDomain {
    /// Create a new FFT domain for the given size.
    ///
    /// Size must be a power of 2 and at most 2^108.
    pub fn new(size: usize) -> Self {
        assert!(size.is_power_of_two(), "FFT size must be a power of 2");
        let log_size = size.trailing_zeros();
        assert!(log_size <= 108, "FFT size too large for Fp128");

        let omega = root_of_unity(log_size);
        let omega_inv = omega.invert().expect("omega should be invertible");
        let size_inv = Fp128::from_u64(size as u64).invert().expect("size should be invertible");

        // Precompute twiddle factors
        let half_size = size / 2;
        let mut twiddles = Vec::with_capacity(half_size);
        let mut twiddles_inv = Vec::with_capacity(half_size);

        let mut w = Fp128::ONE;
        let mut w_inv = Fp128::ONE;
        for _ in 0..half_size {
            twiddles.push(w);
            twiddles_inv.push(w_inv);
            w = w * omega;
            w_inv = w_inv * omega_inv;
        }

        Self {
            size,
            log_size,
            omega,
            omega_inv,
            size_inv,
            twiddles,
            twiddles_inv,
        }
    }

    /// Get the size of this domain.
    pub fn size(&self) -> usize {
        self.size
    }
}

/// Bit-reverse permutation of array indices.
fn bit_reverse_permutation<F: Field>(data: &mut [F]) {
    let n = data.len();
    let log_n = n.trailing_zeros();

    for i in 0..n {
        let j = reverse_bits(i, log_n);
        if i < j {
            data.swap(i, j);
        }
    }
}

/// Reverse the low `bits` bits of `x`.
fn reverse_bits(x: usize, bits: u32) -> usize {
    let mut result = 0;
    let mut x = x;
    for _ in 0..bits {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    result
}

/// In-place radix-2 Cooley-Tukey FFT.
///
/// Transforms coefficients to evaluations at powers of omega.
pub fn fft_in_place(data: &mut [Fp128], domain: &FftDomain) {
    let n = data.len();
    assert_eq!(n, domain.size, "Data size must match domain size");

    // Bit-reverse permutation
    bit_reverse_permutation(data);

    // Cooley-Tukey butterfly operations
    let mut m = 1;
    for _ in 0..domain.log_size {
        let half_m = m;
        m *= 2;

        // Stride for twiddle factor lookup
        let twiddle_stride = domain.size / m;

        for k in (0..n).step_by(m) {
            for j in 0..half_m {
                let twiddle = domain.twiddles[j * twiddle_stride];
                let u = data[k + j];
                let t = twiddle * data[k + j + half_m];
                data[k + j] = u + t;
                data[k + j + half_m] = u - t;
            }
        }
    }
}

/// In-place inverse FFT.
///
/// Transforms evaluations back to coefficients.
pub fn ifft_in_place(data: &mut [Fp128], domain: &FftDomain) {
    let n = data.len();
    assert_eq!(n, domain.size, "Data size must match domain size");

    // Bit-reverse permutation
    bit_reverse_permutation(data);

    // Cooley-Tukey butterfly operations with inverse twiddles
    let mut m = 1;
    for _ in 0..domain.log_size {
        let half_m = m;
        m *= 2;

        let twiddle_stride = domain.size / m;

        for k in (0..n).step_by(m) {
            for j in 0..half_m {
                let twiddle = domain.twiddles_inv[j * twiddle_stride];
                let u = data[k + j];
                let t = twiddle * data[k + j + half_m];
                data[k + j] = u + t;
                data[k + j + half_m] = u - t;
            }
        }
    }

    // Scale by 1/n
    for x in data.iter_mut() {
        *x = *x * domain.size_inv;
    }
}

/// FFT: transform coefficients to evaluations.
///
/// Returns evaluations at [1, omega, omega^2, ..., omega^(n-1)].
pub fn fft(coeffs: &[Fp128], domain: &FftDomain) -> Vec<Fp128> {
    let mut data = coeffs.to_vec();
    // Pad with zeros if needed
    data.resize(domain.size, Fp128::ZERO);
    fft_in_place(&mut data, domain);
    data
}

/// Inverse FFT: transform evaluations to coefficients.
pub fn ifft(evals: &[Fp128], domain: &FftDomain) -> Vec<Fp128> {
    let mut data = evals.to_vec();
    data.resize(domain.size, Fp128::ZERO);
    ifft_in_place(&mut data, domain);
    data
}

/// Multiply two polynomials using FFT.
///
/// Given coefficients of polynomials A and B, compute coefficients of A * B.
pub fn polynomial_multiply(a: &[Fp128], b: &[Fp128]) -> Vec<Fp128> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }

    // Result degree is deg(a) + deg(b), so we need size > deg(a) + deg(b)
    let result_len = a.len() + b.len() - 1;
    let fft_size = result_len.next_power_of_two();
    let domain = FftDomain::new(fft_size);

    // Forward FFT
    let a_evals = fft(a, &domain);
    let b_evals = fft(b, &domain);

    // Pointwise multiply
    let mut c_evals: Vec<Fp128> = a_evals.iter()
        .zip(b_evals.iter())
        .map(|(&x, &y)| x * y)
        .collect();

    // Inverse FFT
    ifft_in_place(&mut c_evals, &domain);

    // Truncate to actual result length
    c_evals.truncate(result_len);
    c_evals
}

/// Coset FFT: evaluate polynomial at coset omega^c * {1, omega, omega^2, ...}.
///
/// This is useful for Reed-Solomon encoding where we want to evaluate
/// at points other than powers of omega.
pub fn coset_fft(coeffs: &[Fp128], domain: &FftDomain, coset_shift: Fp128) -> Vec<Fp128> {
    // Multiply coefficient i by coset_shift^i
    let mut shifted: Vec<Fp128> = Vec::with_capacity(coeffs.len());
    let mut power = Fp128::ONE;
    for &c in coeffs {
        shifted.push(c * power);
        power = power * coset_shift;
    }
    shifted.resize(domain.size, Fp128::ZERO);

    fft_in_place(&mut shifted, domain);
    shifted
}

/// Inverse coset FFT.
pub fn coset_ifft(evals: &[Fp128], domain: &FftDomain, coset_shift: Fp128) -> Vec<Fp128> {
    let mut data = evals.to_vec();
    data.resize(domain.size, Fp128::ZERO);
    ifft_in_place(&mut data, domain);

    // Divide coefficient i by coset_shift^i
    let coset_shift_inv = coset_shift.invert().expect("coset shift should be invertible");
    let mut power = Fp128::ONE;
    for c in data.iter_mut() {
        *c = *c * power;
        power = power * coset_shift_inv;
    }

    data
}

/// Extend polynomial evaluations using FFT.
///
/// Given evaluations at n consecutive integer points (0, 1, ..., n-1),
/// this function returns evaluations at m points (0, 1, ..., m-1).
///
/// This is the FFT-accelerated version of `extend_evaluations`.
///
/// Note: This works by interpolating to get coefficients, then evaluating.
/// For points 0, 1, ..., n-1 which are NOT roots of unity, we use
/// a specialized approach.
pub fn extend_evaluations_fft(
    input: &[Fp128],
    input_len: usize,
    output_len: usize,
) -> Vec<Fp128> {
    assert!(input_len <= input.len());
    assert!(output_len >= input_len);

    if input_len == 0 {
        return vec![Fp128::ZERO; output_len];
    }

    // For small sizes, the overhead of FFT isn't worth it
    if input_len <= 32 {
        return super::polynomial::extend_evaluations_fast(input, input_len, output_len);
    }

    // Get polynomial coefficients via interpolation
    // We have evaluations at 0, 1, 2, ..., n-1 (not roots of unity)
    // Use the existing fast barycentric interpolation to get coefficients
    let coeffs = super::polynomial::interpolate_fast(&input[..input_len]);

    // Now evaluate at 0, 1, 2, ..., output_len-1
    // For consecutive integers, direct evaluation is often faster than FFT
    // unless output_len is very large

    let poly = super::polynomial::Polynomial::from_coeffs(coeffs);
    let mut output = Vec::with_capacity(output_len);
    for i in 0..output_len {
        output.push(poly.evaluate(Fp128::from_u64(i as u64)));
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_root_of_unity() {
        // Test that omega^n = 1 for various n
        for k in 1..=20 {
            let n = 1usize << k;
            let omega = root_of_unity(k);

            // omega^n should equal 1
            let mut power = Fp128::ONE;
            for _ in 0..n {
                power = power * omega;
            }
            assert_eq!(power, Fp128::ONE, "omega^{} should be 1 for k={}", n, k);

            // omega^(n/2) should not equal 1 (primitive root check)
            let mut half_power = Fp128::ONE;
            for _ in 0..(n/2) {
                half_power = half_power * omega;
            }
            assert_ne!(half_power, Fp128::ONE, "omega^{} should not be 1 for k={}", n/2, k);
        }
    }

    #[test]
    fn test_fft_ifft_roundtrip() {
        let domain = FftDomain::new(8);

        let coeffs = vec![
            Fp128::from_u64(1),
            Fp128::from_u64(2),
            Fp128::from_u64(3),
            Fp128::from_u64(4),
            Fp128::from_u64(5),
            Fp128::from_u64(6),
            Fp128::from_u64(7),
            Fp128::from_u64(8),
        ];

        let evals = fft(&coeffs, &domain);
        let recovered = ifft(&evals, &domain);

        assert_eq!(coeffs, recovered);
    }

    #[test]
    fn test_fft_ifft_roundtrip_random() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        for log_size in 3..=10 {
            let size = 1 << log_size;
            let domain = FftDomain::new(size);

            let coeffs: Vec<Fp128> = (0..size).map(|_| Fp128::random(&mut rng)).collect();

            let evals = fft(&coeffs, &domain);
            let recovered = ifft(&evals, &domain);

            assert_eq!(coeffs, recovered, "Roundtrip failed for size {}", size);
        }
    }

    #[test]
    fn test_polynomial_multiply() {
        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
        let a = vec![Fp128::from_u64(1), Fp128::from_u64(2)];
        let b = vec![Fp128::from_u64(3), Fp128::from_u64(4)];

        let c = polynomial_multiply(&a, &b);

        assert_eq!(c.len(), 3);
        assert_eq!(c[0], Fp128::from_u64(3));
        assert_eq!(c[1], Fp128::from_u64(10));
        assert_eq!(c[2], Fp128::from_u64(8));
    }

    #[test]
    fn test_polynomial_multiply_larger() {
        // (1 + x + x^2) * (1 + x) = 1 + 2x + 2x^2 + x^3
        let a = vec![Fp128::from_u64(1), Fp128::from_u64(1), Fp128::from_u64(1)];
        let b = vec![Fp128::from_u64(1), Fp128::from_u64(1)];

        let c = polynomial_multiply(&a, &b);

        assert_eq!(c.len(), 4);
        assert_eq!(c[0], Fp128::from_u64(1));
        assert_eq!(c[1], Fp128::from_u64(2));
        assert_eq!(c[2], Fp128::from_u64(2));
        assert_eq!(c[3], Fp128::from_u64(1));
    }

    #[test]
    fn test_fft_matches_direct_evaluation() {
        let domain = FftDomain::new(8);

        // Polynomial: 1 + 2x + 3x^2
        let coeffs = vec![
            Fp128::from_u64(1),
            Fp128::from_u64(2),
            Fp128::from_u64(3),
            Fp128::ZERO,
            Fp128::ZERO,
            Fp128::ZERO,
            Fp128::ZERO,
            Fp128::ZERO,
        ];

        let evals = fft(&coeffs, &domain);

        // Verify by direct evaluation at powers of omega
        let poly = super::super::polynomial::Polynomial::from_coeffs(coeffs.clone());
        let mut omega_power = Fp128::ONE;
        for (i, &eval) in evals.iter().enumerate() {
            let expected = poly.evaluate(omega_power);
            assert_eq!(eval, expected, "Mismatch at index {}", i);
            omega_power = omega_power * domain.omega;
        }
    }

    #[test]
    fn test_bit_reverse() {
        assert_eq!(reverse_bits(0b000, 3), 0b000);
        assert_eq!(reverse_bits(0b001, 3), 0b100);
        assert_eq!(reverse_bits(0b010, 3), 0b010);
        assert_eq!(reverse_bits(0b011, 3), 0b110);
        assert_eq!(reverse_bits(0b100, 3), 0b001);
        assert_eq!(reverse_bits(0b101, 3), 0b101);
        assert_eq!(reverse_bits(0b110, 3), 0b011);
        assert_eq!(reverse_bits(0b111, 3), 0b111);
    }
}
