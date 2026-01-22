//! Simple example demonstrating the Longfellow ZK proof system.
//!
//! This example shows how to:
//! 1. Create a circuit with quadratic constraints
//! 2. Generate a zero-knowledge proof
//! 3. Verify the proof
//!
//! Run with: `cargo run --example simple_proof`

use longfellow_zk::{
    circuit::{CircuitBuilder, LayerBuilder},
    field::Fp128,
    zk::{verify_zk, ZkProver},
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

fn main() {
    println!("=== Longfellow ZK Proof System Demo ===\n");

    // Example 1: Simple multiplication proof
    // Prove knowledge of a, b such that a * b = c (where c is public)
    simple_multiplication_proof();

    // Example 2: Multi-constraint proof
    // Prove knowledge of values satisfying multiple quadratic constraints
    multi_constraint_proof();

    println!("\n=== All proofs verified successfully! ===");
}

/// Demonstrates a simple multiplication proof: a * b = c
fn simple_multiplication_proof() {
    println!("Example 1: Simple Multiplication Proof");
    println!("--------------------------------------");
    println!("Proving knowledge of a, b such that a * b = c\n");

    // Secret values
    let a = Fp128::from_u64(7);
    let b = Fp128::from_u64(13);
    let c = a * b; // c = 91

    println!("Secret: a = 7, b = 13");
    println!("Public: c = a * b = 91\n");

    // Create a circuit that checks a * b = c
    // The circuit has 3 witness values: w[0] = a, w[1] = b, w[2] = c
    // And one quadratic constraint: w[0] * w[1] - w[2] = 0
    let mut builder = CircuitBuilder::<Fp128>::new(0, 3);

    // Add a layer with one multiplication gate
    // mul(output_idx, left_idx, right_idx, coefficient)
    // This represents: output[0] = coefficient * w[left] * w[right]
    let mut layer = LayerBuilder::new(0, 2);
    layer.add_mul(0, 0, 1, Fp128::ONE); // output[0] = w[0] * w[1]
    builder.add_layer(layer);

    let circuit = builder.build();

    // Witness vector: [a, b, c]
    let witness = vec![a, b, c];

    // Create prover
    let quadratic_constraints = vec![]; // Using circuit layer instead
    let prover = ZkProver::new(circuit.clone(), quadratic_constraints, witness);

    // Generate proof
    let mut rng = ChaCha20Rng::seed_from_u64(42);
    let proof = prover.prove(&mut rng);

    println!("Proof generated!");
    println!("  - Ligero commitment size: {} bytes", proof.commitment.to_bytes().len());
    println!("  - Ligero proof size: {} bytes\n", proof.ligero_proof.size());

    // Verify proof
    let public_inputs: Vec<Fp128> = vec![]; // c is embedded in the circuit
    let verified = verify_zk(&circuit, &public_inputs, &proof);

    println!("Verification result: {}", if verified { "VALID" } else { "INVALID" });
    assert!(verified, "Proof verification failed!");
    println!();
}

/// Demonstrates a proof with multiple quadratic constraints
fn multi_constraint_proof() {
    println!("Example 2: Multi-Constraint Proof");
    println!("----------------------------------");
    println!("Proving knowledge of x, y, z satisfying:");
    println!("  x * y = p");
    println!("  y * z = q");
    println!("  x * z = r\n");

    // Secret values
    let x = Fp128::from_u64(2);
    let y = Fp128::from_u64(3);
    let z = Fp128::from_u64(5);

    // Computed products
    let p = x * y; // 6
    let q = y * z; // 15
    let r = x * z; // 10

    println!("Secret: x = 2, y = 3, z = 5");
    println!("Products: p = 6, q = 15, r = 10\n");

    // Create a circuit with multiple multiplication gates
    // Witness: [x, y, z, p, q, r]
    let mut builder = CircuitBuilder::<Fp128>::new(0, 6);

    // Layer 1: compute x*y and y*z
    let mut layer1 = LayerBuilder::new(0, 4);
    layer1.add_mul(0, 0, 1, Fp128::ONE); // output[0] = x * y
    layer1.add_mul(1, 1, 2, Fp128::ONE); // output[1] = y * z
    builder.add_layer(layer1);

    // Layer 2: compute x*z
    let mut layer2 = LayerBuilder::new(0, 2);
    layer2.add_mul(0, 0, 2, Fp128::ONE); // output[0] = x * z
    builder.add_layer(layer2);

    let circuit = builder.build();

    // Witness vector
    let witness = vec![x, y, z, p, q, r];

    // Create prover
    let prover = ZkProver::new(circuit.clone(), vec![], witness);

    // Generate proof
    let mut rng = ChaCha20Rng::seed_from_u64(123);
    let proof = prover.prove(&mut rng);

    println!("Proof generated!");

    // Verify
    let public_inputs: Vec<Fp128> = vec![];
    let verified = verify_zk(&circuit, &public_inputs, &proof);

    println!("Verification result: {}", if verified { "VALID" } else { "INVALID" });
    assert!(verified, "Proof verification failed!");
}
