//! Ligero prover implementation.
//!
//! This module implements the Ligero prover that commits to a witness
//! and generates proofs for linear and quadratic constraints.

use crate::field::Field;
use crate::merkle::{hash_data, MerkleTree};
use crate::polynomial::{axpy, extend, pointwise_mul, pointwise_sub};
use crate::transcript::Transcript;
use rand::Rng;

use super::commitment::LigeroCommitment;
use super::params::LigeroParams;
use super::proof::LigeroProof;

/// A quadratic constraint (x, y, z) meaning W[x] * W[y] = W[z].
#[derive(Clone, Copy, Debug)]
pub struct QuadraticConstraint {
    /// Index of the first multiplicand in the witness vector.
    pub x: usize,
    /// Index of the second multiplicand in the witness vector.
    pub y: usize,
    /// Index of the product in the witness vector.
    pub z: usize,
}

/// A linear constraint term (witness_idx, constraint_idx, coefficient).
#[derive(Clone, Copy, Debug)]
pub struct LinearConstraintTerm {
    /// Index of the witness element this term refers to.
    pub witness_idx: usize,
    /// Index of the constraint this term belongs to.
    pub constraint_idx: usize,
    /// Coefficient value for this term.
    pub coefficient: usize,
}

/// The Ligero prover.
pub struct LigeroProver<F: Field> {
    /// Parameters for the commitment scheme.
    params: LigeroParams,

    /// The witness vector.
    witness: Vec<F>,

    /// The tableau matrix (NROW x NCOL).
    tableau: Vec<Vec<F>>,

    /// Quadratic constraints.
    quadratic_constraints: Vec<QuadraticConstraint>,
}

impl<F: Field> LigeroProver<F> {
    /// Create a new Ligero prover for the given witness and constraints.
    pub fn new(
        witness: Vec<F>,
        quadratic_constraints: Vec<QuadraticConstraint>,
        security_level: usize,
    ) -> Self {
        let params =
            LigeroParams::new(witness.len(), quadratic_constraints.len(), security_level);

        Self {
            params,
            witness,
            tableau: Vec::new(),
            quadratic_constraints,
        }
    }

    /// Create a prover with explicit parameters (for test vectors).
    pub fn with_params(
        witness: Vec<F>,
        quadratic_constraints: Vec<QuadraticConstraint>,
        params: LigeroParams,
    ) -> Self {
        Self {
            params,
            witness,
            tableau: Vec::new(),
            quadratic_constraints,
        }
    }

    /// Commit to the witness vector.
    ///
    /// Returns the commitment and prepares internal state for proving.
    pub fn commit<R: Rng>(&mut self, rng: &mut R) -> LigeroCommitment {
        // Initialize tableau
        self.tableau = vec![vec![F::ZERO; self.params.ncol]; self.params.nrow];

        // Layout the tableau
        self.layout_zk_rows(rng);
        self.layout_witness_rows(rng);
        self.layout_quadratic_rows(rng);

        // Build Merkle tree from columns
        let num_leaves = self.params.ncol - self.params.dblock;
        let mut tree = MerkleTree::new(num_leaves);

        for j in self.params.dblock..self.params.ncol {
            // Hash the column
            let mut col_bytes = Vec::new();
            for i in 0..self.params.nrow {
                col_bytes.extend_from_slice(&self.tableau[i][j].to_bytes());
            }
            let leaf = hash_data(&col_bytes);
            tree.set_leaf(j - self.params.dblock, leaf);
        }

        let root = tree.build();
        LigeroCommitment::new(root)
    }

    /// Generate a proof for the given constraints.
    pub fn prove<R: Rng>(
        &self,
        transcript: &mut Transcript,
        commitment: &LigeroCommitment,
        linear_constraints: &[LinearConstraintTerm],
        b: &[F],
        _rng: &mut R,
    ) -> LigeroProof<F> {
        let params = &self.params;

        // Write commitment to transcript
        transcript.write_bytes(&commitment.to_bytes());

        // Generate challenge for low-degree test
        let u: Vec<F> = transcript.generate_challenges(params.nrow);

        // Compute LDT response: sum of u[i] * T[i][0..BLOCK]
        let mut ldt = vec![F::ZERO; params.block];
        for i in 0..params.nrow {
            for j in 0..params.block {
                ldt[j] = ldt[j] + u[i] * self.tableau[i][j];
            }
        }

        // Write LDT to transcript
        for elem in &ldt {
            transcript.write_field_element(elem);
        }

        // Generate challenges for linear and quadratic constraints
        let alpha_l: Vec<F> = transcript.generate_challenges(b.len().max(1));
        let alpha_q: Vec<Vec<F>> = (0..params.nq)
            .map(|_| transcript.generate_challenges(3))
            .collect();

        // Compute inner product vector A
        let a = self.inner_product_vector(linear_constraints, &alpha_l, &alpha_q);

        // Compute dot proof
        let dot = self.dot_proof(&a);

        // Write dot to transcript
        for elem in &dot {
            transcript.write_field_element(elem);
        }

        // Generate quadratic challenge
        let u_quad: Vec<F> = transcript.generate_challenges(params.nqt);

        // Compute quadratic proof
        let qpr = self.quadratic_proof(&u_quad);

        // Write quadratic proof to transcript
        for elem in &qpr {
            transcript.write_field_element(elem);
        }

        // Generate challenge indices
        let challenge_indices = transcript
            .generate_indices_without_replacement(params.ncol - params.dblock, params.nreq);

        // Collect requested columns
        let columns = self.requested_columns(&challenge_indices);

        // Build Merkle proof
        let num_leaves = params.ncol - params.dblock;
        let mut tree = MerkleTree::new(num_leaves);

        for j in params.dblock..params.ncol {
            let mut col_bytes = Vec::new();
            for i in 0..params.nrow {
                col_bytes.extend_from_slice(&self.tableau[i][j].to_bytes());
            }
            let leaf = hash_data(&col_bytes);
            tree.set_leaf(j - params.dblock, leaf);
        }
        tree.build();

        let merkle_proof = tree.compressed_proof(&challenge_indices);

        LigeroProof::new(ldt, dot, qpr, columns, merkle_proof)
    }

    /// Get the parameters.
    pub fn params(&self) -> &LigeroParams {
        &self.params
    }

    /// Layout the zero-knowledge rows (ILDT, IDOT, IQD).
    fn layout_zk_rows<R: Rng>(&mut self, rng: &mut R) {
        let params = &self.params;

        // Row 0 (ILDT): random row
        let random_row: Vec<F> = (0..params.block).map(|_| F::random(rng)).collect();
        self.tableau[0] = extend(&random_row, params.block, params.ncol);

        // Row 1 (IDOT): random row with sum constraint
        // Use block-size to ensure all rows have the same polynomial degree
        let mut z: Vec<F> = (0..params.block).map(|_| F::random(rng)).collect();
        // Ensure sum of z[NREQ..NREQ+WR] = 0
        let mut sum = F::ZERO;
        let end_idx = (params.nreq + params.wr).min(params.block);
        for i in params.nreq..end_idx.saturating_sub(1) {
            sum = sum + z[i];
        }
        if end_idx > 0 && end_idx - 1 < params.block {
            z[end_idx - 1] = -sum;
        }
        self.tableau[1] = extend(&z, params.block, params.ncol);

        // Row 2 (IQD): random row with witness positions zeroed
        // Use block-size to ensure all rows have the same polynomial degree
        let mut zq: Vec<F> = (0..params.block).map(|_| F::random(rng)).collect();
        let end_idx = (params.nreq + params.wr).min(params.block);
        for i in params.nreq..end_idx {
            zq[i] = F::ZERO;
        }
        self.tableau[2] = extend(&zq, params.block, params.ncol);
    }

    /// Layout the witness rows.
    fn layout_witness_rows<R: Rng>(&mut self, rng: &mut R) {
        let params = &self.params;

        for i in 0..params.nwrow {
            let mut row = Vec::with_capacity(params.block);

            // Random padding
            for _ in 0..params.nreq {
                row.push(F::random(rng));
            }

            // Witness values
            for j in 0..params.wr {
                let w_idx = i * params.wr + j;
                if w_idx < self.witness.len() {
                    row.push(self.witness[w_idx]);
                } else {
                    row.push(F::ZERO);
                }
            }

            // Pad to BLOCK size if needed
            while row.len() < params.block {
                row.push(F::ZERO);
            }

            self.tableau[params.iw + i] = extend(&row, params.block, params.ncol);
        }
    }

    /// Layout the quadratic constraint rows.
    fn layout_quadratic_rows<R: Rng>(&mut self, rng: &mut R) {
        let params = &self.params;

        for i in 0..params.nqt {
            let mut qx = Vec::with_capacity(params.block);
            let mut qy = Vec::with_capacity(params.block);
            let mut qz = Vec::with_capacity(params.block);

            // Random padding
            for _ in 0..params.nreq {
                qx.push(F::random(rng));
                qy.push(F::random(rng));
                qz.push(F::random(rng));
            }

            // Quadratic constraint values
            for j in 0..params.wr {
                let q_idx = i * params.qr + j;
                if q_idx < self.quadratic_constraints.len() {
                    let qc = &self.quadratic_constraints[q_idx];
                    let wx = self.witness.get(qc.x).copied().unwrap_or(F::ZERO);
                    let wy = self.witness.get(qc.y).copied().unwrap_or(F::ZERO);
                    let wz = self.witness.get(qc.z).copied().unwrap_or(F::ZERO);

                    qx.push(wx);
                    qy.push(wy);
                    qz.push(wz);
                } else {
                    qx.push(F::ZERO);
                    qy.push(F::ZERO);
                    qz.push(F::ZERO);
                }
            }

            // Pad to BLOCK size
            while qx.len() < params.block {
                qx.push(F::ZERO);
                qy.push(F::ZERO);
                qz.push(F::ZERO);
            }

            let row_base = params.iq + 3 * i;
            self.tableau[row_base] = extend(&qx, params.block, params.ncol);
            self.tableau[row_base + 1] = extend(&qy, params.block, params.ncol);
            self.tableau[row_base + 2] = extend(&qz, params.block, params.ncol);
        }
    }

    /// Compute the inner product vector for linear and quadratic constraints.
    fn inner_product_vector(
        &self,
        linear: &[LinearConstraintTerm],
        alpha_l: &[F],
        alpha_q: &[Vec<F>],
    ) -> Vec<F> {
        let params = &self.params;
        let total_size = params.nqw * params.wr;
        let mut a = vec![F::ZERO; total_size];

        // Process linear constraints
        for term in linear {
            if term.witness_idx < self.witness.len() && term.constraint_idx < alpha_l.len() {
                let coeff = F::from_u64(term.coefficient as u64);
                a[term.witness_idx] =
                    a[term.witness_idx] + alpha_l[term.constraint_idx] * coeff;
            }
        }

        // Process quadratic constraints
        let a_x_base = params.nwrow * params.wr;
        let a_y_base = a_x_base + params.nqt * params.wr;
        let a_z_base = a_y_base + params.nqt * params.wr;

        for (i, qc) in self.quadratic_constraints.iter().enumerate() {
            if i >= alpha_q.len() {
                break;
            }

            let triple_idx = i / params.qr;
            let offset = i % params.qr;
            let ia = triple_idx * params.wr + offset;

            // Copy constraints: Qx[ia] - W[x] = 0
            if a_x_base + ia < a.len() {
                a[a_x_base + ia] = a[a_x_base + ia] + alpha_q[i][0];
            }
            if qc.x < a.len() {
                a[qc.x] = a[qc.x] - alpha_q[i][0];
            }

            if a_y_base + ia < a.len() {
                a[a_y_base + ia] = a[a_y_base + ia] + alpha_q[i][1];
            }
            if qc.y < a.len() {
                a[qc.y] = a[qc.y] - alpha_q[i][1];
            }

            if a_z_base + ia < a.len() {
                a[a_z_base + ia] = a[a_z_base + ia] + alpha_q[i][2];
            }
            if qc.z < a.len() {
                a[qc.z] = a[qc.z] - alpha_q[i][2];
            }
        }

        a
    }

    /// Compute the dot product proof.
    fn dot_proof(&self, a: &[F]) -> Vec<F> {
        let params = &self.params;
        let mut y = self.tableau[params.idot()][..params.dblock].to_vec();

        for i in 0..params.nqw {
            let mut a_ext = vec![F::ZERO; params.block];
            for j in 0..params.wr {
                let idx = i * params.wr + j;
                if idx < a.len() {
                    a_ext[params.nreq + j] = a[idx];
                }
            }

            let af = extend(&a_ext, params.block, params.dblock);

            // y += Af .* T[i + IW][0..DBLOCK]
            let row_idx = params.iw + i;
            if row_idx < self.tableau.len() {
                for j in 0..params.dblock {
                    y[j] = y[j] + af[j] * self.tableau[row_idx][j];
                }
            }
        }

        y
    }

    /// Compute the quadratic proof.
    fn quadratic_proof(&self, u_quad: &[F]) -> Vec<F> {
        let params = &self.params;
        let mut y = self.tableau[params.iqd()][..params.dblock].to_vec();

        for i in 0..params.nqt {
            let iqx = params.iq + 3 * i;
            let iqy = iqx + 1;
            let iqz = iqx + 2;

            if iqz >= self.tableau.len() {
                break;
            }

            // tmp = z[i] - x[i] * y[i]
            let z_row = &self.tableau[iqz][..params.dblock];
            let x_row = &self.tableau[iqx][..params.dblock];
            let y_row = &self.tableau[iqy][..params.dblock];

            let xy = pointwise_mul(x_row, y_row);
            let tmp = pointwise_sub(z_row, &xy);

            // y += u_quad[i] * tmp
            if i < u_quad.len() {
                y = axpy(&y, u_quad[i], &tmp);
            }
        }

        y
    }

    /// Collect requested columns.
    fn requested_columns(&self, indices: &[usize]) -> Vec<Vec<F>> {
        let params = &self.params;

        indices
            .iter()
            .map(|&i| {
                let col_idx = i + params.dblock;
                (0..params.nrow)
                    .map(|row| self.tableau[row][col_idx])
                    .collect()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp128;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_commitment() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let witness: Vec<Fp128> = vec![
            Fp128::from_u64(1),
            Fp128::from_u64(2),
            Fp128::from_u64(3),
            Fp128::from_u64(4),
        ];

        // 1 * 2 = 2 (stored at witness[4] if it existed, but we simplify)
        let constraints = vec![];

        let mut prover = LigeroProver::new(witness, constraints, 128);
        let commitment = prover.commit(&mut rng);

        // Commitment should be a 32-byte hash
        assert_eq!(commitment.to_bytes().len(), 32);
    }
}
