//! Multi-field example demonstrating proofs over the different fields
//! supported by the Longfellow ZK proof system.
//!
//! This example shows the same circuit proven over three different fields:
//! - Fp128: prime field 2^128 - 2^108 + 1 (FieldID 0x06)
//! - GF(2^128): binary extension field (FieldID 0x04)
//! - Fp256: NIST P-256 base field (FieldID 0x01)
//!
//! The GF(2^128) and Fp256 fields are required by the IETF spec for the
//! full Longfellow anonymous credential scheme (ECDSA signature proofs).
//!
//! Run with: `cargo run --example multi_field`

use longfellow_zk::{
    circuit::{CircuitBuilder, LayerBuilder},
    field::{Field, FieldId, Fp128, Fp256, GF2_128},
    zk::{verify_zk, ZkProver},
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

fn main() {
    println!("=== Longfellow ZK Multi-Field Demo ===\n");

    println!("Demonstrating proofs over all fields required by IETF spec:\n");

    demo_field::<Fp128>("Fp128 (prime 2^128 - 2^108 + 1)", 7, 13);
    demo_field::<GF2_128>("GF(2^128) (binary extension field)", 7, 13);
    demo_field::<Fp256>("Fp256 (NIST P-256 base field)", 7, 13);

    println!("\n--- Binary Field Characteristics ---\n");
    demo_binary_field_properties();

    println!("\n--- P-256 Field Properties ---\n");
    demo_p256_properties();

    println!("\n=== All multi-field demos completed! ===");
}

/// Generic demo that works over any field implementing the Field trait.
fn demo_field<F: Field>(name: &str, a: u64, b: u64) {
    println!("Field: {}", name);
    println!("  FieldID: {:#04x}", F::FIELD_ID as u8);
    println!("  Byte length: {}", F::BYTE_LEN);
    println!("  Characteristic: {}", if F::CHARACTERISTIC == 0 {
        "large prime".to_string()
    } else {
        F::CHARACTERISTIC.to_string()
    });

    // Create a simple multiplication circuit: c = a * b
    let mut builder = CircuitBuilder::<F>::new(0, 3);
    let mut layer = LayerBuilder::new(0, 2);
    layer.add_mul(0, 0, 1, F::ONE);
    builder.add_layer(layer);
    let circuit = builder.build();

    // Compute witness: a * b = c (in the field!)
    let fa = F::from_u64(a);
    let fb = F::from_u64(b);
    let fc = fa * fb;

    println!("  Proving: {} * {} = {:?}", a, b, fc);
    println!("  (Note: in GF(2^128), multiplication is carryless polynomial mult)");

    let witness = vec![fa, fb, fc];
    let prover = ZkProver::new(circuit.clone(), vec![], witness);

    let mut rng = ChaCha20Rng::seed_from_u64(42);
    let proof = prover.prove(&mut rng);

    let verified = verify_zk(&circuit, &[], &proof);
    println!("  Proof verified: {}", if verified { "YES" } else { "NO" });
    println!("  Proof size: {} bytes\n", proof.size());
    assert!(verified);
}

/// Demonstrate the unique properties of the binary extension field.
fn demo_binary_field_properties() {
    // In GF(2^128), characteristic is 2, so 1 + 1 = 0
    let one = GF2_128::ONE;
    println!("  1 + 1 = {:?} (characteristic 2)", one + one);
    assert_eq!(one + one, GF2_128::ZERO);

    // Negation is the identity: -a = a
    let a = GF2_128::from_u64(42);
    println!("  -42 = {:?} (negation is identity)", -a);
    assert_eq!(-a, a);

    // The P2 constant for sumcheck is X, not 2
    println!("  P2 (sumcheck eval point) = X = {:?}", <GF2_128 as Field>::P2);
    assert_eq!(<GF2_128 as Field>::P2, GF2_128::X);
    // Note: this is NOT the same as "two" since 1+1=0 in this field!
    assert_ne!(<GF2_128 as Field>::P2, GF2_128::ONE + GF2_128::ONE);

    // x^128 reduces to x^7 + x^2 + x + 1 = 0x87
    let x = GF2_128::X;
    let x128 = x.pow(128);
    println!("  x^128 mod Q(x) = {:#x} (= x^7 + x^2 + x + 1)", x128.to_u128());
    assert_eq!(x128, GF2_128::from_u64(0x87));

    // Multiplication is polynomial multiplication mod Q(x), NOT integer mult
    // 7 * 13 in integers = 91 = 0x5b
    // 7 * 13 in GF(2^128) = (x^2+x+1)(x^3+x^2+1)
    //   = x^5 + x^4 + x^2  (from x^2 term)
    //   + x^4 + x^3 + x    (from x term)
    //   + x^3 + x^2 + 1    (from 1 term)
    //   = x^5 + x + 1 = 0b100011 = 35 = 0x23
    //   (coefficients cancel pairwise since 1+1=0)
    let result = GF2_128::from_u64(7) * GF2_128::from_u64(13);
    println!("  7 * 13 = {:#x} (carryless: differs from integer 91)", result.to_u128());
    assert_eq!(result, GF2_128::from_u64(0b100011));

    // The full proof system works correctly over this field
    println!("  GF(2^128) enables efficient SHA-256 circuits in Longfellow ZK");
}

/// Demonstrate NIST P-256 base field properties.
fn demo_p256_properties() {
    // The P-256 prime: 2^256 - 2^224 + 2^192 + 2^96 - 1
    println!("  Prime p = 2^256 - 2^224 + 2^192 + 2^96 - 1");
    println!("  This is the base field of the NIST P-256 elliptic curve");
    println!("  Used in Longfellow for proving ECDSA-P256 signatures");

    // Verify: p - 1 + 1 = 0
    let p_minus_1 = Fp256::ZERO - Fp256::ONE;
    assert_eq!(p_minus_1 + Fp256::ONE, Fp256::ZERO);
    println!("  (p - 1) + 1 = 0: verified");

    // Inversion works via Fermat's little theorem
    let a = Fp256::from_u64(12345);
    let a_inv = a.invert().unwrap();
    assert_eq!(a * a_inv, Fp256::ONE);
    println!("  12345 * 12345^(-1) = 1: verified");

    println!("  FieldID: {:#04x} (matches spec value for p256)", FieldId::P256 as u8);
}
