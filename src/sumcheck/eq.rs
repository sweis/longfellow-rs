//! EQ array and binding operations for the sumcheck protocol.
//!
//! The EQ_{n}[i, j] array is defined as 1 if i = j and i < n, 0 otherwise.
//! These operations are used for computing bindings in the sumcheck protocol.

use crate::field::Field;

/// Bind a single variable in an array.
///
/// Given array A and field element x, returns B such that:
/// B[i] = (1 - x) * A[2 * i] + x * A[2 * i + 1]
pub fn bind<F: Field>(a: &[F], x: F) -> Vec<F> {
    let n = (a.len() + 1) / 2;
    let mut b = Vec::with_capacity(n);

    let one_minus_x = F::ONE - x;

    for i in 0..n {
        let a0 = a.get(2 * i).copied().unwrap_or(F::ZERO);
        let a1 = a.get(2 * i + 1).copied().unwrap_or(F::ZERO);
        b.push(one_minus_x * a0 + x * a1);
    }

    b
}

/// Bind multiple variables in sequence.
///
/// bindv(A, X) = A if X is empty, otherwise bindv(bind(A, X[0]), X[1..])
pub fn bindv<F: Field>(a: &[F], x: &[F]) -> Vec<F> {
    let mut result = a.to_vec();
    for &xi in x {
        result = bind(&result, xi);
    }
    result
}

/// Compute bindv(EQ_n, X) efficiently.
///
/// For n = 2^l and X of size l, this computes the binding of the EQ array
/// without explicitly constructing the full EQ array.
pub fn bind_eq<F: Field>(l: usize, x: &[F]) -> Vec<F> {
    assert_eq!(x.len(), l, "X must have exactly l elements");

    if l == 0 {
        return vec![F::ONE];
    }

    let n = 1 << l;
    let mut b = vec![F::ZERO; n];

    // Recursive construction
    let a = bind_eq(l - 1, &x[1..]);

    let one_minus_x0 = F::ONE - x[0];

    for i in 0..a.len() {
        b[2 * i] = one_minus_x0 * a[i];
        b[2 * i + 1] = x[0] * a[i];
    }

    b
}

/// Compute bindv(EQ_m, X) for arbitrary m (not necessarily a power of 2).
///
/// Pads m to the next power of 2 and ignores extra elements.
#[allow(dead_code)]
pub fn bind_eq_general<F: Field>(m: usize, x: &[F]) -> Vec<F> {
    if m == 0 {
        return vec![];
    }

    // Compute l such that 2^l >= m
    let l = (usize::BITS - m.leading_zeros()) as usize;

    if l != x.len() {
        // Need to handle the case where x is the wrong size
        // For now, we'll pad or truncate
        let padded_x: Vec<F> = x
            .iter()
            .copied()
            .chain(std::iter::repeat(F::ZERO))
            .take(l)
            .collect();

        let full = bind_eq(l, &padded_x);
        return full[..m].to_vec();
    }

    let full = bind_eq(l, x);
    full[..m].to_vec()
}

/// Bind the first dimension of a 2D array.
#[allow(dead_code)]
pub fn bind_2d<F: Field>(a: &[Vec<F>], x: F) -> Vec<Vec<F>> {
    let n = (a.len() + 1) / 2;
    let mut b = Vec::with_capacity(n);

    let one_minus_x = F::ONE - x;

    for i in 0..n {
        let a0 = a.get(2 * i);
        let a1 = a.get(2 * i + 1);

        let width = a0.map(|v| v.len()).unwrap_or(0).max(a1.map(|v| v.len()).unwrap_or(0));
        let mut row = Vec::with_capacity(width);

        for j in 0..width {
            let v0 = a0.and_then(|v| v.get(j)).copied().unwrap_or(F::ZERO);
            let v1 = a1.and_then(|v| v.get(j)).copied().unwrap_or(F::ZERO);
            row.push(one_minus_x * v0 + x * v1);
        }

        b.push(row);
    }

    b
}

/// Bind the first dimension of a 3D array (QUAD).
#[allow(dead_code)]
pub fn bind_3d<F: Field>(q: &[Vec<Vec<F>>], x: F) -> Vec<Vec<Vec<F>>> {
    let n = (q.len() + 1) / 2;
    let mut b = Vec::with_capacity(n);

    let one_minus_x = F::ONE - x;

    for i in 0..n {
        let q0 = q.get(2 * i);
        let q1 = q.get(2 * i + 1);

        let l_size = q0.map(|v| v.len()).unwrap_or(0).max(q1.map(|v| v.len()).unwrap_or(0));
        let mut plane = Vec::with_capacity(l_size);

        for l in 0..l_size {
            let row0 = q0.and_then(|v| v.get(l));
            let row1 = q1.and_then(|v| v.get(l));

            let r_size = row0.map(|v| v.len()).unwrap_or(0).max(row1.map(|v| v.len()).unwrap_or(0));
            let mut row = Vec::with_capacity(r_size);

            for r in 0..r_size {
                let v0 = row0.and_then(|v| v.get(r)).copied().unwrap_or(F::ZERO);
                let v1 = row1.and_then(|v| v.get(r)).copied().unwrap_or(F::ZERO);
                row.push(one_minus_x * v0 + x * v1);
            }

            plane.push(row);
        }

        b.push(plane);
    }

    b
}

/// Transpose a 2D array.
pub fn transpose<F: Field>(a: &[Vec<F>]) -> Vec<Vec<F>> {
    if a.is_empty() {
        return vec![];
    }

    let rows = a.len();
    let cols = a.iter().map(|r| r.len()).max().unwrap_or(0);

    let mut result = vec![vec![F::ZERO; rows]; cols];

    for (i, row) in a.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            result[j][i] = val;
        }
    }

    result
}

/// Transpose a 2D QUAD array (swap l and r dimensions).
pub fn transpose_quad<F: Field>(q: &[Vec<F>]) -> Vec<Vec<F>> {
    transpose(q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp128;

    #[test]
    fn test_bind() {
        let a = vec![
            Fp128::from_u64(1),
            Fp128::from_u64(2),
            Fp128::from_u64(3),
            Fp128::from_u64(4),
        ];
        let x = Fp128::from_u64(0);

        let b = bind(&a, x);
        // When x = 0, B[i] = A[2*i]
        assert_eq!(b[0], Fp128::from_u64(1));
        assert_eq!(b[1], Fp128::from_u64(3));

        let x = Fp128::from_u64(1);
        let b = bind(&a, x);
        // When x = 1, B[i] = A[2*i + 1]
        assert_eq!(b[0], Fp128::from_u64(2));
        assert_eq!(b[1], Fp128::from_u64(4));
    }

    #[test]
    fn test_bind_eq() {
        // EQ_4 is a 4x4 identity-like array
        // bindv(EQ_4, [0, 0]) should give [1, 0, 0, 0]
        let result = bind_eq(2, &[Fp128::ZERO, Fp128::ZERO]);
        assert_eq!(result[0], Fp128::ONE);
        assert_eq!(result[1], Fp128::ZERO);
        assert_eq!(result[2], Fp128::ZERO);
        assert_eq!(result[3], Fp128::ZERO);

        // bindv(EQ_4, [1, 0]) should give [0, 1, 0, 0]
        let result = bind_eq(2, &[Fp128::ONE, Fp128::ZERO]);
        assert_eq!(result[0], Fp128::ZERO);
        assert_eq!(result[1], Fp128::ONE);
        assert_eq!(result[2], Fp128::ZERO);
        assert_eq!(result[3], Fp128::ZERO);
    }

    #[test]
    fn test_transpose() {
        let a = vec![
            vec![Fp128::from_u64(1), Fp128::from_u64(2)],
            vec![Fp128::from_u64(3), Fp128::from_u64(4)],
        ];

        let t = transpose(&a);
        assert_eq!(t[0][0], Fp128::from_u64(1));
        assert_eq!(t[0][1], Fp128::from_u64(3));
        assert_eq!(t[1][0], Fp128::from_u64(2));
        assert_eq!(t[1][1], Fp128::from_u64(4));
    }
}
