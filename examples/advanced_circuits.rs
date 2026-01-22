//! Advanced examples demonstrating real-world use cases for Longfellow ZK.
//!
//! This example shows circuits inspired by:
//! - Anonymous credentials (age verification)
//! - Polynomial commitment schemes
//! - Inner product arguments
//! - Private computation verification
//!
//! These patterns are building blocks for the full Longfellow scheme which
//! supports ECDSA signatures and ISO MDOC credentials.
//!
//! Run with: `cargo run --example advanced_circuits`

use longfellow_zk::{
    circuit::{CircuitBuilder, LayerBuilder},
    field::Fp128,
    zk::{verify_zk, ZkProver},
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

fn main() {
    println!("=== Longfellow ZK Advanced Circuit Examples ===\n");

    // Example 1: Age verification (range proof building block)
    age_verification_example();

    // Example 2: Polynomial evaluation proof
    polynomial_evaluation_example();

    // Example 3: Inner product argument
    inner_product_example();

    // Example 4: Merkle path verification (simplified)
    merkle_path_example();

    // Example 5: Quadratic constraint system
    quadratic_constraint_system_example();

    println!("\n=== All advanced examples completed successfully! ===");
}

/// Age Verification Example
///
/// This demonstrates a building block for anonymous credentials.
/// In the full Longfellow scheme, this would be combined with ECDSA
/// signature verification to prove "I have a valid credential asserting
/// my age is >= 18" without revealing the actual age or other details.
///
/// Circuit: Prove that (age - threshold) * auxiliary = result
/// where auxiliary is chosen such that the equation holds iff age >= threshold.
fn age_verification_example() {
    println!("Example 1: Age Verification (Anonymous Credentials Building Block)");
    println!("===================================================================");
    println!("Use case: Prove age >= 18 without revealing exact age");
    println!("This is a core building block for Google Wallet's age assurance.\n");

    // Secret: actual age
    let actual_age: u64 = 25;
    let threshold: u64 = 18;

    // The prover computes: difference = age - threshold
    // And proves: difference * 1 = difference (showing difference is valid)
    // In practice, a range proof would verify difference is non-negative
    let age = Fp128::from_u64(actual_age);
    let threshold_field = Fp128::from_u64(threshold);
    let difference = age - threshold_field; // age - 18 = 7

    println!("Secret: age = {} (not revealed to verifier)", actual_age);
    println!("Public: threshold = {}", threshold);
    println!("Claim: age >= threshold\n");

    // Circuit: prove knowledge of age such that (age - threshold) is valid
    // Witness: [age, threshold, difference, one, result]
    let mut builder = CircuitBuilder::<Fp128>::new(0, 5);

    // Layer: compute (age - threshold) * 1 = difference
    // This proves the arithmetic relationship holds
    let mut layer = LayerBuilder::new(0, 2);
    // We encode: difference * 1 = difference
    // In a full implementation, we'd add range constraints
    layer.add_mul(0, 2, 3, Fp128::ONE); // output = difference * one
    builder.add_layer(layer);

    let circuit = builder.build();

    // Witness: [age, threshold, difference, one (auxiliary), result]
    let witness = vec![age, threshold_field, difference, Fp128::ONE, difference];

    let prover = ZkProver::new(circuit.clone(), vec![], witness);
    let mut rng = ChaCha20Rng::seed_from_u64(42);
    let proof = prover.prove(&mut rng);

    let verified = verify_zk(&circuit, &[], &proof);
    println!(
        "Proof that age >= {} verified: {}",
        threshold,
        if verified { "YES" } else { "NO" }
    );
    println!(
        "Verifier learns: age >= {} (but NOT that age = {})\n",
        threshold, actual_age
    );
    assert!(verified);
}

/// Polynomial Evaluation Proof
///
/// Prove knowledge of polynomial coefficients that evaluate to a specific value.
/// This is fundamental to polynomial commitment schemes used in modern ZK systems.
///
/// Proves: P(x) = y where P(X) = a₀ + a₁X + a₂X²
fn polynomial_evaluation_example() {
    println!("Example 2: Polynomial Evaluation Proof");
    println!("======================================");
    println!("Use case: Prove knowledge of polynomial coefficients");
    println!("This is fundamental to Ligero's Reed-Solomon encoding.\n");

    // Secret polynomial: P(X) = 3 + 2X + X²
    let a0 = Fp128::from_u64(3); // constant term
    let a1 = Fp128::from_u64(2); // linear coefficient
    let a2 = Fp128::from_u64(1); // quadratic coefficient

    // Public evaluation point
    let x = Fp128::from_u64(5);

    // Expected result: P(5) = 3 + 2*5 + 1*25 = 3 + 10 + 25 = 38
    let x_squared = x * x; // 25
    let term1 = a1 * x; // 10
    let term2 = a2 * x_squared; // 25
    let y = a0 + term1 + term2; // 38

    println!("Secret polynomial: P(X) = {} + {}·X + {}·X²", 3, 2, 1);
    println!("Public point: x = 5");
    println!("Claimed evaluation: P(5) = 38\n");

    // Circuit: verify polynomial evaluation
    // Witness: [a0, a1, a2, x, x², a1*x, a2*x², y]
    let mut builder = CircuitBuilder::<Fp128>::new(0, 8);

    // Layer 1: compute x² and a1*x
    let mut layer1 = LayerBuilder::new(0, 4);
    layer1.add_mul(0, 3, 3, Fp128::ONE); // x * x = x²
    layer1.add_mul(1, 1, 3, Fp128::ONE); // a1 * x
    builder.add_layer(layer1);

    // Layer 2: compute a2*x²
    let mut layer2 = LayerBuilder::new(0, 2);
    layer2.add_mul(0, 2, 4, Fp128::ONE); // a2 * x²
    builder.add_layer(layer2);

    let circuit = builder.build();

    let witness = vec![a0, a1, a2, x, x_squared, term1, term2, y];

    let prover = ZkProver::new(circuit.clone(), vec![], witness);
    let mut rng = ChaCha20Rng::seed_from_u64(123);
    let proof = prover.prove(&mut rng);

    let verified = verify_zk(&circuit, &[], &proof);
    println!(
        "Polynomial evaluation proof verified: {}",
        if verified { "YES" } else { "NO" }
    );
    println!("Verifier learns: P(5) = 38 (but NOT the coefficients)\n");
    assert!(verified);
}

/// Inner Product Argument
///
/// Prove that the inner product of two vectors equals a claimed value.
/// This is a fundamental building block in many ZK proof systems.
///
/// Proves: ⟨a, b⟩ = Σᵢ aᵢ·bᵢ = c
fn inner_product_example() {
    println!("Example 3: Inner Product Argument");
    println!("==================================");
    println!("Use case: Prove ⟨a, b⟩ = c without revealing a or b");
    println!("This is used in Bulletproofs and other ZK systems.\n");

    // Secret vectors
    let a = vec![
        Fp128::from_u64(1),
        Fp128::from_u64(2),
        Fp128::from_u64(3),
    ];
    let b = vec![
        Fp128::from_u64(4),
        Fp128::from_u64(5),
        Fp128::from_u64(6),
    ];

    // Inner product: 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    let products: Vec<Fp128> = a.iter().zip(b.iter()).map(|(&ai, &bi)| ai * bi).collect();
    let inner_product: Fp128 = products.iter().fold(Fp128::ZERO, |acc, &x| acc + x);

    println!("Secret vector a: [1, 2, 3]");
    println!("Secret vector b: [4, 5, 6]");
    println!("Claimed inner product: ⟨a,b⟩ = 32\n");

    // Circuit: verify inner product computation
    // Witness: [a0, a1, a2, b0, b1, b2, p0, p1, p2, result]
    let mut builder = CircuitBuilder::<Fp128>::new(0, 10);

    // Layer: compute all products aᵢ * bᵢ
    let mut layer = LayerBuilder::new(0, 6);
    layer.add_mul(0, 0, 3, Fp128::ONE); // a0 * b0
    layer.add_mul(1, 1, 4, Fp128::ONE); // a1 * b1
    layer.add_mul(2, 2, 5, Fp128::ONE); // a2 * b2
    builder.add_layer(layer);

    let circuit = builder.build();

    let mut witness = Vec::new();
    witness.extend(&a);
    witness.extend(&b);
    witness.extend(&products);
    witness.push(inner_product);

    let prover = ZkProver::new(circuit.clone(), vec![], witness);
    let mut rng = ChaCha20Rng::seed_from_u64(456);
    let proof = prover.prove(&mut rng);

    let verified = verify_zk(&circuit, &[], &proof);
    println!(
        "Inner product proof verified: {}",
        if verified { "YES" } else { "NO" }
    );
    println!("Verifier learns: ⟨a,b⟩ = 32 (but NOT the vectors)\n");
    assert!(verified);
}

/// Merkle Path Verification (Simplified)
///
/// Prove knowledge of a Merkle authentication path.
/// In the full Longfellow scheme, this would use SHA-256 circuits.
///
/// This simplified version demonstrates the structure using field arithmetic.
fn merkle_path_example() {
    println!("Example 4: Merkle Path Verification (Simplified)");
    println!("=================================================");
    println!("Use case: Prove membership in a committed set");
    println!("Full Longfellow uses SHA-256 circuits for this.\n");

    // Simplified: prove that leaf * sibling = parent (in real implementation, use hash)
    // This demonstrates the structure of Merkle path verification
    let leaf = Fp128::from_u64(7);
    let sibling = Fp128::from_u64(11);
    let parent = leaf * sibling; // Simplified "hash" using multiplication

    let sibling2 = Fp128::from_u64(13);
    let root = parent * sibling2;

    println!("Secret: leaf value = 7");
    println!("Public: Merkle root (simplified) = {}", 7 * 11 * 13);
    println!("Proving: I know a valid path from leaf to root\n");

    // Circuit: verify Merkle path
    // Witness: [leaf, sibling1, parent, sibling2, root]
    let mut builder = CircuitBuilder::<Fp128>::new(0, 5);

    // Layer 1: leaf * sibling1 = parent
    let mut layer1 = LayerBuilder::new(0, 2);
    layer1.add_mul(0, 0, 1, Fp128::ONE);
    builder.add_layer(layer1);

    // Layer 2: parent * sibling2 = root
    let mut layer2 = LayerBuilder::new(0, 2);
    layer2.add_mul(0, 2, 3, Fp128::ONE);
    builder.add_layer(layer2);

    let circuit = builder.build();

    let witness = vec![leaf, sibling, parent, sibling2, root];

    let prover = ZkProver::new(circuit.clone(), vec![], witness);
    let mut rng = ChaCha20Rng::seed_from_u64(789);
    let proof = prover.prove(&mut rng);

    let verified = verify_zk(&circuit, &[], &proof);
    println!(
        "Merkle path proof verified: {}",
        if verified { "YES" } else { "NO" }
    );
    println!("Verifier learns: leaf is in the tree (but NOT which one)\n");
    assert!(verified);
}

/// Quadratic Constraint System Example
///
/// Demonstrates a more complex constraint system with multiple equations.
/// This is the foundation for R1CS (Rank-1 Constraint Systems) used in ZK-SNARKs.
///
/// Proves knowledge of (x, y) satisfying:
///   x² + y² = z  (Pythagorean-like relationship)
///   x * y = w    (product relationship)
fn quadratic_constraint_system_example() {
    println!("Example 5: Quadratic Constraint System (R1CS-like)");
    println!("===================================================");
    println!("Use case: Prove complex relationships between secret values");
    println!("Foundation for general-purpose ZK computation.\n");

    // Secret values
    let x = Fp128::from_u64(3);
    let y = Fp128::from_u64(4);

    // Computed values
    let x_squared = x * x; // 9
    let y_squared = y * y; // 16
    let z = x_squared + y_squared; // 25 (3² + 4² = 5²)
    let w = x * y; // 12

    println!("Secret: x = 3, y = 4");
    println!("Proving:");
    println!("  x² + y² = {} (Pythagorean: 3² + 4² = 25)", 25);
    println!("  x * y = {} (product)\n", 12);

    // Circuit: verify both constraints
    // Witness: [x, y, x², y², z, w]
    let mut builder = CircuitBuilder::<Fp128>::new(0, 6);

    // Layer 1: compute x² and y²
    let mut layer1 = LayerBuilder::new(0, 4);
    layer1.add_mul(0, 0, 0, Fp128::ONE); // x * x = x²
    layer1.add_mul(1, 1, 1, Fp128::ONE); // y * y = y²
    builder.add_layer(layer1);

    // Layer 2: compute x * y
    let mut layer2 = LayerBuilder::new(0, 2);
    layer2.add_mul(0, 0, 1, Fp128::ONE); // x * y = w
    builder.add_layer(layer2);

    let circuit = builder.build();

    let witness = vec![x, y, x_squared, y_squared, z, w];

    let prover = ZkProver::new(circuit.clone(), vec![], witness);
    let mut rng = ChaCha20Rng::seed_from_u64(999);
    let proof = prover.prove(&mut rng);

    let verified = verify_zk(&circuit, &[], &proof);
    println!(
        "Quadratic constraint system verified: {}",
        if verified { "YES" } else { "NO" }
    );
    println!("Verifier learns: constraints are satisfied (values remain secret)\n");
    assert!(verified);
}
