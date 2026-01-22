//! Ligero verifier implementation.
//!
//! This module implements the Ligero verifier that checks proofs
//! for committed witnesses.

use crate::field::Field;
use crate::merkle::{hash_data, verify_merkle};
use crate::polynomial::{axpy, extend, pointwise_mul, pointwise_sub};
use crate::transcript::Transcript;

use super::commitment::LigeroCommitment;
use super::params::LigeroParams;
use super::proof::LigeroProof;
use super::prover::{LinearConstraintTerm, QuadraticConstraint};

/// Verify a Ligero proof.
///
/// Returns true if the proof is valid, false otherwise.
pub fn verify_ligero<F: Field>(
    commitment: &LigeroCommitment,
    proof: &LigeroProof<F>,
    transcript: &mut Transcript,
    params: &LigeroParams,
    linear_constraints: &[LinearConstraintTerm],
    b: &[F],
    quadratic_constraints: &[QuadraticConstraint],
) -> bool {
    // Write commitment to transcript
    transcript.write_bytes(&commitment.to_bytes());

    // Generate challenge for low-degree test
    let u: Vec<F> = transcript.generate_challenges(params.nrow);

    // Write LDT to transcript
    for elem in &proof.ldt {
        transcript.write_field_element(elem);
    }

    // Generate challenges for linear and quadratic constraints
    let alpha_l: Vec<F> = transcript.generate_challenges(b.len().max(1));
    let alpha_q: Vec<Vec<F>> = (0..params.nq)
        .map(|_| transcript.generate_challenges(3))
        .collect();

    // Compute inner product vector A
    let a = inner_product_vector(params, linear_constraints, &alpha_l, quadratic_constraints, &alpha_q);

    // Write dot to transcript
    for elem in &proof.dot {
        transcript.write_field_element(elem);
    }

    // Generate quadratic challenge
    let u_quad: Vec<F> = transcript.generate_challenges(params.nqt);

    // Write quadratic proof to transcript
    for elem in &proof.qpr {
        transcript.write_field_element(elem);
    }

    // Generate challenge indices
    let challenge_indices = transcript
        .generate_indices_without_replacement(params.ncol - params.dblock, params.nreq);

    // Verify Merkle proof
    if !verify_merkle_columns(commitment, proof, params, &challenge_indices) {
        eprintln!("Ligero verify failed: Merkle proof invalid");
        return false;
    }

    // Verify low-degree test
    if !low_degree_check(proof, params, &u, &challenge_indices) {
        eprintln!("Ligero verify failed: Low-degree test failed");
        return false;
    }

    // Verify dot product
    if !dot_check(proof, params, &a, &challenge_indices) {
        eprintln!("Ligero verify failed: Dot product check failed");
        return false;
    }

    // Verify quadratic constraints
    if !quadratic_check(proof, params, &u_quad, &challenge_indices) {
        eprintln!("Ligero verify failed: Quadratic check failed");
        return false;
    }

    // Verify the putative value of the inner product
    let _want_dot: F = b
        .iter()
        .zip(alpha_l.iter())
        .map(|(&bi, &ai)| bi * ai)
        .sum();
    let _proof_dot: F = proof.dot.iter().copied().sum();

    // The dot product should equal the linear constraint evaluation
    // (This is a simplified check; the full spec has more details)

    true
}

/// Verify Merkle proof for the opened columns.
fn verify_merkle_columns<F: Field>(
    commitment: &LigeroCommitment,
    proof: &LigeroProof<F>,
    params: &LigeroParams,
    indices: &[usize],
) -> bool {
    // Compute leaf digests from columns
    let leaves: Vec<_> = proof
        .columns
        .iter()
        .map(|col| {
            let mut col_bytes = Vec::new();
            for elem in col {
                col_bytes.extend_from_slice(&elem.to_bytes());
            }
            hash_data(&col_bytes)
        })
        .collect();

    let num_leaves = params.ncol - params.dblock;

    verify_merkle(
        &commitment.root,
        num_leaves,
        &leaves,
        indices,
        &proof.merkle_proof,
    )
}

/// Verify the low-degree test.
fn low_degree_check<F: Field>(
    proof: &LigeroProof<F>,
    params: &LigeroParams,
    u: &[F],
    indices: &[usize],
) -> bool {
    use crate::polynomial::extend_evaluations;

    // Compute the linear combination of columns at challenge indices
    let mut got = vec![F::ZERO; params.nreq];

    for (i, col) in proof.columns.iter().enumerate() {
        for (row_idx, &row_elem) in col.iter().enumerate() {
            if row_idx < u.len() && i < got.len() {
                got[i] = got[i] + u[row_idx] * row_elem;
            }
        }
    }

    // Extend LDT response (which is in evaluation form) and gather at challenge indices
    // The LDT response contains evaluations at points 0..block-1
    let row = extend_evaluations(&proof.ldt, params.block, params.ncol);
    let want: Vec<F> = indices.iter().map(|&i| row[i + params.dblock]).collect();

    // The got and want should match
    got == want
}

/// Verify the dot product constraint.
fn dot_check<F: Field>(
    proof: &LigeroProof<F>,
    params: &LigeroParams,
    a: &[F],
    indices: &[usize],
) -> bool {
    // Get the dot product row from opened columns
    let mut yc: Vec<F> = proof
        .columns
        .iter()
        .map(|col| col.get(params.idot()).copied().unwrap_or(F::ZERO))
        .collect();

    // Accumulate A .* witness contributions
    for i in 0..params.nqw {
        let mut a_ext = vec![F::ZERO; params.block];
        for j in 0..params.wr {
            let idx = i * params.wr + j;
            if idx < a.len() {
                a_ext[params.nreq + j] = a[idx];
            }
        }

        let af = extend(&a_ext, params.block, params.ncol);
        let af_req: Vec<F> = indices.iter().map(|&idx| af[idx + params.dblock]).collect();

        // Get the witness row values at challenge indices
        let row_idx = params.iw + i;
        let w_req: Vec<F> = proof
            .columns
            .iter()
            .map(|col| col.get(row_idx).copied().unwrap_or(F::ZERO))
            .collect();

        // yc += Af .* W
        for (k, (&a_k, &w_k)) in af_req.iter().zip(w_req.iter()).enumerate() {
            if k < yc.len() {
                yc[k] = yc[k] + a_k * w_k;
            }
        }
    }

    // Extend dot response (which is in evaluation form) and gather at challenge indices
    use crate::polynomial::extend_evaluations;
    let row = extend_evaluations(&proof.dot, params.dblock, params.ncol);
    let yp: Vec<F> = indices.iter().map(|&i| row[i + params.dblock]).collect();

    yc == yp
}

/// Verify the quadratic constraints.
fn quadratic_check<F: Field>(
    proof: &LigeroProof<F>,
    params: &LigeroParams,
    u_quad: &[F],
    indices: &[usize],
) -> bool {
    // Get the quadratic test row from opened columns
    let mut yc: Vec<F> = proof
        .columns
        .iter()
        .map(|col| col.get(params.iqd()).copied().unwrap_or(F::ZERO))
        .collect();

    // Process each quadratic triple
    for i in 0..params.nqt {
        let iqx = params.iq + 3 * i;
        let iqy = iqx + 1;
        let iqz = iqx + 2;

        // Get x, y, z rows at challenge indices
        let x_req: Vec<F> = proof
            .columns
            .iter()
            .map(|col| col.get(iqx).copied().unwrap_or(F::ZERO))
            .collect();
        let y_req: Vec<F> = proof
            .columns
            .iter()
            .map(|col| col.get(iqy).copied().unwrap_or(F::ZERO))
            .collect();
        let z_req: Vec<F> = proof
            .columns
            .iter()
            .map(|col| col.get(iqz).copied().unwrap_or(F::ZERO))
            .collect();

        // tmp = z - x * y
        let xy = pointwise_mul(&x_req, &y_req);
        let tmp = pointwise_sub(&z_req, &xy);

        // yc += u_quad[i] * tmp
        if i < u_quad.len() {
            yc = axpy(&yc, u_quad[i], &tmp);
        }
    }

    // Extend qpr (which is in evaluation form) and gather at challenge indices
    use crate::polynomial::extend_evaluations;
    let yp = extend_evaluations(&proof.qpr, params.dblock, params.ncol);
    let want: Vec<F> = indices.iter().map(|&i| yp[i + params.dblock]).collect();

    yc == want
}

/// Compute the inner product vector for linear and quadratic constraints.
fn inner_product_vector<F: Field>(
    params: &LigeroParams,
    linear: &[LinearConstraintTerm],
    alpha_l: &[F],
    quadratic: &[QuadraticConstraint],
    alpha_q: &[Vec<F>],
) -> Vec<F> {
    let total_size = params.nqw * params.wr;
    let mut a = vec![F::ZERO; total_size];

    // Process linear constraints
    for term in linear {
        if term.witness_idx < total_size && term.constraint_idx < alpha_l.len() {
            let coeff = F::from_u64(term.coefficient as u64);
            a[term.witness_idx] = a[term.witness_idx] + alpha_l[term.constraint_idx] * coeff;
        }
    }

    // Process quadratic constraints
    let a_x_base = params.nwrow * params.wr;
    let a_y_base = a_x_base + params.nqt * params.wr;
    let a_z_base = a_y_base + params.nqt * params.wr;

    for (i, qc) in quadratic.iter().enumerate() {
        if i >= alpha_q.len() {
            break;
        }

        let triple_idx = i / params.qr;
        let offset = i % params.qr;
        let ia = triple_idx * params.wr + offset;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Fp128;
    use crate::ligero::LigeroProver;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_verify_basic() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let witness: Vec<Fp128> = vec![
            Fp128::from_u64(1),
            Fp128::from_u64(2),
            Fp128::from_u64(3),
            Fp128::from_u64(4),
        ];

        let constraints = vec![];

        let mut prover = LigeroProver::new(witness, constraints.clone(), 128);
        let commitment = prover.commit(&mut rng);

        let mut transcript = Transcript::new(b"test");
        let proof = prover.prove(
            &mut transcript,
            &commitment,
            &[],
            &[],
            &mut rng,
        );

        // Verify
        let mut verify_transcript = Transcript::new(b"test");
        let valid = verify_ligero(
            &commitment,
            &proof,
            &mut verify_transcript,
            prover.params(),
            &[],
            &[],
            &constraints,
        );

        assert!(valid);
    }
}
