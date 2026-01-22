//! Polynomial operations for the Longfellow ZK scheme.
//!
//! This module implements polynomial interpolation, evaluation, and extension
//! (Reed-Solomon encoding) used in the Ligero commitment scheme.

use crate::field::Field;

/// A polynomial represented by its coefficients.
///
/// Coefficients are stored in ascending order of degree:
/// poly[0] + poly[1]*x + poly[2]*x^2 + ...
#[derive(Clone, Debug)]
pub struct Polynomial<F: Field> {
    /// Coefficients in ascending order of degree.
    pub coeffs: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    /// Create a polynomial from coefficients.
    pub fn from_coeffs(coeffs: Vec<F>) -> Self {
        Self { coeffs }
    }

    /// Create the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: vec![] }
    }

    /// Get the degree of the polynomial (-1 for zero polynomial).
    pub fn degree(&self) -> isize {
        if self.coeffs.is_empty() {
            -1
        } else {
            (self.coeffs.len() - 1) as isize
        }
    }

    /// Evaluate the polynomial at a point.
    pub fn evaluate(&self, x: F) -> F {
        if self.coeffs.is_empty() {
            return F::ZERO;
        }

        // Horner's method
        let mut result = self.coeffs[self.coeffs.len() - 1];
        for i in (0..self.coeffs.len() - 1).rev() {
            result = result * x + self.coeffs[i];
        }
        result
    }

    /// Interpolate a polynomial from points and values.
    ///
    /// Given points x[0], ..., x[n-1] and values y[0], ..., y[n-1],
    /// returns the unique polynomial of degree < n such that P(x[i]) = y[i].
    pub fn interpolate(points: &[F], values: &[F]) -> Self {
        assert_eq!(points.len(), values.len());
        let n = points.len();

        if n == 0 {
            return Self::zero();
        }

        // Lagrange interpolation
        let mut result = vec![F::ZERO; n];

        for i in 0..n {
            // Compute the i-th Lagrange basis polynomial
            let mut basis = vec![F::ZERO; n];
            basis[0] = F::ONE;
            let mut degree = 0;

            for j in 0..n {
                if i == j {
                    continue;
                }

                // Multiply by (x - points[j]) / (points[i] - points[j])
                let denom = (points[i] - points[j]).invert().unwrap();
                let factor = -points[j] * denom;

                // Shift and add
                for k in (1..=degree + 1).rev() {
                    basis[k] = basis[k] * factor + basis[k - 1] * denom;
                }
                basis[0] = basis[0] * factor;
                degree += 1;
            }

            // Add values[i] * basis to result
            for (k, &b) in basis.iter().enumerate() {
                result[k] = result[k] + values[i] * b;
            }
        }

        // Remove trailing zeros
        while result.last() == Some(&F::ZERO) && result.len() > 1 {
            result.pop();
        }

        Self { coeffs: result }
    }
}

/// Extend an array using polynomial evaluation (Reed-Solomon encoding).
///
/// Given `input` of length `input_len`, this function treats the input as
/// polynomial coefficients and evaluates the polynomial at points
/// `0, 1, 2, ..., output_len - 1` to produce the output.
///
/// This is the core operation for Ligero's low-degree test.
pub fn extend<F: Field>(input: &[F], input_len: usize, output_len: usize) -> Vec<F> {
    assert!(input_len <= input.len());
    assert!(output_len >= input_len);

    // Treat input[0..input_len] as polynomial coefficients
    let poly = Polynomial::from_coeffs(input[..input_len].to_vec());

    // Evaluate at 0, 1, 2, ..., output_len - 1
    let mut output = Vec::with_capacity(output_len);
    for i in 0..output_len {
        let x = F::from_u64(i as u64);
        output.push(poly.evaluate(x));
    }

    output
}

/// Extend an array treating input as evaluations at consecutive points.
///
/// Given evaluations at points 0, 1, ..., input_len - 1, interpolate the
/// polynomial and evaluate at points input_len, ..., output_len - 1.
pub fn extend_evaluations<F: Field>(
    input: &[F],
    input_len: usize,
    output_len: usize,
) -> Vec<F> {
    assert!(input_len <= input.len());
    assert!(output_len >= input_len);

    if input_len == 0 {
        return vec![F::ZERO; output_len];
    }

    // Construct evaluation points
    let points: Vec<F> = (0..input_len).map(|i| F::from_u64(i as u64)).collect();
    let values = &input[..input_len];

    // Interpolate
    let poly = Polynomial::interpolate(&points, values);

    // Evaluate at all output points
    let mut output = Vec::with_capacity(output_len);
    for i in 0..output_len {
        let x = F::from_u64(i as u64);
        output.push(poly.evaluate(x));
    }

    output
}

/// Multiply two polynomials (represented as evaluation arrays).
///
/// Given arrays representing polynomial evaluations, compute the
/// element-wise (Hadamard) product.
pub fn pointwise_mul<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect()
}

/// Add two polynomial evaluation arrays.
pub fn pointwise_add<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
}

/// Subtract polynomial evaluation arrays.
pub fn pointwise_sub<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(&x, &y)| x - y).collect()
}

/// Scale a polynomial evaluation array.
pub fn pointwise_scale<F: Field>(a: &[F], scalar: F) -> Vec<F> {
    a.iter().map(|&x| x * scalar).collect()
}

/// Compute axpy: result = a + scalar * b
pub fn axpy<F: Field>(a: &[F], scalar: F, b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| x + scalar * y)
        .collect()
}

/// Gather elements from an array at the specified indices.
pub fn gather<F: Field>(array: &[F], indices: &[usize]) -> Vec<F> {
    indices.iter().map(|&i| array[i]).collect()
}

/// Compute the sum of an array.
pub fn sum<F: Field>(array: &[F]) -> F {
    array.iter().copied().sum()
}

/// Compute the inner product of two arrays.
pub fn inner_product<F: Field>(a: &[F], b: &[F]) -> F {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

/// Lagrange interpolation coefficient at a single point.
///
/// Computes the value of the i-th Lagrange basis polynomial at point x.
pub fn lagrange_basis<F: Field>(points: &[F], i: usize, x: F) -> F {
    let n = points.len();
    let mut result = F::ONE;

    for j in 0..n {
        if i != j {
            let denom = points[i] - points[j];
            let num = x - points[j];
            result = result * num * denom.invert().unwrap();
        }
    }

    result
}

/// Compute Lagrange coefficients for interpolation at three points P0=0, P1=1, P2.
///
/// Returns the coefficients [lag_0(x), lag_1(x), lag_2(x)] such that
/// P(x) = lag_0(x) * P(0) + lag_1(x) * P(1) + lag_2(x) * P(2)
pub fn lagrange_coefficients_3<F: Field>(x: F, p2: F) -> [F; 3] {
    // lag_0(x) = (x - 1)(x - P2) / ((0 - 1)(0 - P2)) = (x - 1)(x - P2) / P2
    // lag_1(x) = (x - 0)(x - P2) / ((1 - 0)(1 - P2)) = x(x - P2) / (1 - P2)
    // lag_2(x) = (x - 0)(x - 1) / ((P2 - 0)(P2 - 1)) = x(x - 1) / (P2(P2 - 1))

    let one = F::ONE;
    let _zero = F::ZERO;

    let p2_inv = p2.invert().unwrap();
    let one_minus_p2_inv = (one - p2).invert().unwrap();
    let p2_p2m1_inv = (p2 * (p2 - one)).invert().unwrap();

    let lag_0 = (x - one) * (x - p2) * p2_inv;
    let lag_1 = x * (x - p2) * one_minus_p2_inv;
    let lag_2 = x * (x - one) * p2_p2m1_inv;

    [lag_0, lag_1, lag_2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp128;

    #[test]
    fn test_polynomial_evaluation() {
        // P(x) = 1 + 2x + 3x^2
        let poly = Polynomial::from_coeffs(vec![
            Fp128::from_u64(1),
            Fp128::from_u64(2),
            Fp128::from_u64(3),
        ]);

        // P(0) = 1
        assert_eq!(poly.evaluate(Fp128::ZERO), Fp128::from_u64(1));

        // P(1) = 1 + 2 + 3 = 6
        assert_eq!(poly.evaluate(Fp128::ONE), Fp128::from_u64(6));

        // P(2) = 1 + 4 + 12 = 17
        assert_eq!(poly.evaluate(Fp128::from_u64(2)), Fp128::from_u64(17));
    }

    #[test]
    fn test_polynomial_interpolation() {
        // Interpolate through (0, 1), (1, 6), (2, 17)
        // Should give P(x) = 1 + 2x + 3x^2
        let points = vec![
            Fp128::from_u64(0),
            Fp128::from_u64(1),
            Fp128::from_u64(2),
        ];
        let values = vec![
            Fp128::from_u64(1),
            Fp128::from_u64(6),
            Fp128::from_u64(17),
        ];

        let poly = Polynomial::interpolate(&points, &values);

        // Verify the interpolation
        for (i, &p) in points.iter().enumerate() {
            assert_eq!(poly.evaluate(p), values[i]);
        }
    }

    #[test]
    fn test_extend() {
        // Extend a linear polynomial [1, 2] (representing 1 + 2x)
        let input = vec![Fp128::from_u64(1), Fp128::from_u64(2)];
        let output = extend(&input, 2, 4);

        // Expected: P(0) = 1, P(1) = 3, P(2) = 5, P(3) = 7
        assert_eq!(output[0], Fp128::from_u64(1));
        assert_eq!(output[1], Fp128::from_u64(3));
        assert_eq!(output[2], Fp128::from_u64(5));
        assert_eq!(output[3], Fp128::from_u64(7));
    }

    #[test]
    fn test_extend_evaluations() {
        // Given evaluations at 0, 1, 2
        let input = vec![
            Fp128::from_u64(1),
            Fp128::from_u64(6),
            Fp128::from_u64(17),
        ];
        let output = extend_evaluations(&input, 3, 5);

        // The polynomial is 1 + 2x + 3x^2
        // P(3) = 1 + 6 + 27 = 34
        // P(4) = 1 + 8 + 48 = 57
        assert_eq!(output[3], Fp128::from_u64(34));
        assert_eq!(output[4], Fp128::from_u64(57));
    }

    #[test]
    fn test_lagrange_basis() {
        let points = vec![
            Fp128::from_u64(0),
            Fp128::from_u64(1),
            Fp128::from_u64(2),
        ];

        // L_0(0) = 1, L_1(0) = 0, L_2(0) = 0
        assert_eq!(lagrange_basis(&points, 0, Fp128::from_u64(0)), Fp128::ONE);
        assert_eq!(lagrange_basis(&points, 1, Fp128::from_u64(0)), Fp128::ZERO);
        assert_eq!(lagrange_basis(&points, 2, Fp128::from_u64(0)), Fp128::ZERO);

        // L_0(1) = 0, L_1(1) = 1, L_2(1) = 0
        assert_eq!(lagrange_basis(&points, 0, Fp128::from_u64(1)), Fp128::ZERO);
        assert_eq!(lagrange_basis(&points, 1, Fp128::from_u64(1)), Fp128::ONE);
        assert_eq!(lagrange_basis(&points, 2, Fp128::from_u64(1)), Fp128::ZERO);
    }
}

#[cfg(test)]
mod tests_large {
    use super::*;
    use crate::field::Fp128;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_interpolation_large() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        // Create a random polynomial of degree < 23
        let n = 23;
        let coeffs: Vec<Fp128> = (0..n).map(|_| Fp128::random(&mut rng)).collect();
        let poly = Polynomial::from_coeffs(coeffs.clone());

        // Evaluate at points 0..n-1
        let points: Vec<Fp128> = (0..n).map(|i| Fp128::from_u64(i as u64)).collect();
        let values: Vec<Fp128> = points.iter().map(|&p| poly.evaluate(p)).collect();

        // Interpolate
        let interp_poly = Polynomial::interpolate(&points, &values);

        // Check that interpolated poly has same evaluations at original points
        for (i, &p) in points.iter().enumerate() {
            let original = values[i];
            let interpolated = interp_poly.evaluate(p);
            assert_eq!(original, interpolated, "Mismatch at point {}", i);
        }

        // Check at a new point (e.g., 74)
        let test_point = Fp128::from_u64(74);
        let original_at_74 = poly.evaluate(test_point);
        let interpolated_at_74 = interp_poly.evaluate(test_point);

        println!("Original poly at 74: {:?}", original_at_74);
        println!("Interpolated poly at 74: {:?}", interpolated_at_74);

        assert_eq!(original_at_74, interpolated_at_74, "Mismatch at point 74");
    }

    #[test]
    fn test_extend_extend_evaluations_consistency() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let block = 23;
        let ncol = 92;

        // Create random coefficients
        let coeffs: Vec<Fp128> = (0..block).map(|_| Fp128::random(&mut rng)).collect();

        // Extend (coefficients -> evaluations)
        let extended = extend(&coeffs, block, ncol);

        // Now use extend_evaluations on first block evaluations
        let re_extended = extend_evaluations(&extended, block, ncol);

        // They should match at all points
        for j in 0..ncol {
            if extended[j] != re_extended[j] {
                println!("Mismatch at j={}: extended={:?}, re_extended={:?}", j, extended[j], re_extended[j]);
                assert_eq!(extended[j], re_extended[j], "Mismatch at point {}", j);
            }
        }
        println!("All {} points match between extend and extend_evaluations", ncol);
    }

    #[test]
    fn test_linear_combination_interpolation() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let block = 23;
        let ncol = 92;
        let nrow = 4;

        // Create nrow random polynomials (as coefficients)
        let poly_coeffs: Vec<Vec<Fp128>> = (0..nrow)
            .map(|_| (0..block).map(|_| Fp128::random(&mut rng)).collect())
            .collect();

        // Extend each polynomial
        let polys: Vec<Vec<Fp128>> = poly_coeffs
            .iter()
            .map(|c| extend(c, block, ncol))
            .collect();

        // Random linear combination coefficients
        let u: Vec<Fp128> = (0..nrow).map(|_| Fp128::random(&mut rng)).collect();

        // Compute linear combination at all points (direct method)
        let mut direct: Vec<Fp128> = vec![Fp128::ZERO; ncol];
        for j in 0..ncol {
            for i in 0..nrow {
                direct[j] = direct[j] + u[i] * polys[i][j];
            }
        }

        // Compute linear combination at first block points only
        let ldt: Vec<Fp128> = (0..block)
            .map(|j| {
                let mut sum = Fp128::ZERO;
                for i in 0..nrow {
                    sum = sum + u[i] * polys[i][j];
                }
                sum
            })
            .collect();

        // Interpolate from ldt and extend
        let interpolated = extend_evaluations(&ldt, block, ncol);

        // Compare at all points
        let mut mismatches = 0;
        for j in 0..ncol {
            if direct[j] != interpolated[j] {
                if mismatches < 5 {
                    println!("Mismatch at j={}: direct={:?}, interpolated={:?}", j, direct[j], interpolated[j]);
                }
                mismatches += 1;
            }
        }
        println!("Total mismatches: {}/{}", mismatches, ncol);
        assert_eq!(mismatches, 0, "Found {} mismatches", mismatches);
    }
}
