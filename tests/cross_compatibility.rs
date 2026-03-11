//! Cross-compatibility tests for the Longfellow ZK implementation.
//!
//! These tests verify that our implementation produces results consistent
//! with the reference implementations:
//! - google/longfellow-zk (C++)
//! - abetterinternet/zk-cred-longfellow (Rust)

use longfellow_zk::field::{Field, FieldId, Fp128, Fp256, GF2_128};
use longfellow_zk::hash::{HashFunction, Sha256Hash};
use longfellow_zk::merkle::{hash_pair, MerkleTreeGeneric};
use longfellow_zk::polynomial::{extend, Polynomial};
use longfellow_zk::transcript::Transcript;
use longfellow_zk::{LigeroParams, LigeroProof};

/// Test that the field modulus is correct.
/// p = 2^128 - 2^108 + 1
#[test]
fn test_field_modulus() {
    // The modulus p = 2^128 - 2^108 + 1
    // In hex: 0xfffff00000000000_0000000000000001
    // p - 1 = 0xfffff00000000000_0000000000000000

    // Test: p - 1 should be the maximum valid element
    // And (p-1) + 1 should wrap to 0
    // p - 1 = [0x0000000000000000, 0xfffff00000000000]
    let max = Fp128::from_raw([0x0000000000000000, 0xfffff00000000000]);
    let one = Fp128::ONE;
    let result = max + one;
    assert_eq!(result, Fp128::ZERO, "p - 1 + 1 should equal 0 mod p");
}

/// Test the key field identity: 2^128 ≡ 2^108 - 1 (mod p)
#[test]
fn test_field_reduction_identity() {
    // 2^128 mod p should equal 2^108 - 1
    // Since 2^128 - 2^108 + 1 ≡ 0 (mod p), we have 2^128 ≡ 2^108 - 1 (mod p)

    // Compute 2^108 - 1 = 0x00000fffffffffff_ffffffffffffffff
    let two_108_minus_1 = Fp128::from_raw([0xffffffffffffffff, 0x00000fffffffffff]);

    // Compute 2^64 as a field element
    let two_64 = Fp128::from_raw([0, 1]);

    // 2^128 = (2^64)^2
    let two_128 = two_64 * two_64;

    assert_eq!(
        two_128, two_108_minus_1,
        "2^128 should equal 2^108 - 1 mod p"
    );
}

/// Test field serialization format (16 bytes, little-endian).
#[test]
fn test_field_serialization() {
    // Test vector: 42 as a field element
    let val = Fp128::from_u64(42);
    let bytes = val.to_bytes();

    // Should be 16 bytes
    assert_eq!(bytes.len(), 16, "Field element should serialize to 16 bytes");

    // In little-endian: 42 = 0x2a in lowest byte
    assert_eq!(bytes[0], 0x2a, "Lowest byte should be 0x2a");
    for i in 1..16 {
        assert_eq!(bytes[i], 0, "Higher bytes should be 0");
    }

    // Verify round-trip
    let restored = Fp128::from_bytes(&bytes);
    assert_eq!(restored, val, "Round-trip serialization should preserve value");
}

/// Test that -1 in the field equals p - 1.
#[test]
fn test_field_negation() {
    let one = Fp128::ONE;
    let neg_one = -one;

    // -1 should equal p - 1 = [0x0000000000000000, 0xfffff00000000000]
    let p_minus_1 = Fp128::from_raw([0x0000000000000000, 0xfffff00000000000]);

    assert_eq!(neg_one, p_minus_1, "-1 should equal p - 1");

    // Also verify: neg_one + one = 0
    assert_eq!(neg_one + one, Fp128::ZERO, "-1 + 1 should equal 0");
}

/// Test field inversion using Fermat's little theorem.
#[test]
fn test_field_inversion() {
    // Test vectors
    let test_values: Vec<u64> = vec![1, 2, 3, 7, 13, 42, 1000, 65537];

    for v in test_values {
        let x = Fp128::from_u64(v);
        let x_inv = x.invert().expect("non-zero element should be invertible");
        let product = x * x_inv;
        assert_eq!(product, Fp128::ONE, "{} * {}^(-1) should equal 1", v, v);
    }
}

/// Test SHA-256 hash consistency for Merkle tree.
#[test]
fn test_sha256_hash_consistency() {
    // Hash of empty data using SHA-256 explicitly
    let empty_hash = Sha256Hash::hash(&[]);
    assert_eq!(
        empty_hash.len(),
        32,
        "SHA-256 output should be 32 bytes"
    );

    // Known SHA-256 test vector
    // SHA256("abc") = ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad
    let abc_hash = Sha256Hash::hash(b"abc");
    let expected: [u8; 32] = [
        0xba, 0x78, 0x16, 0xbf, 0x8f, 0x01, 0xcf, 0xea,
        0x41, 0x41, 0x40, 0xde, 0x5d, 0xae, 0x22, 0x23,
        0xb0, 0x03, 0x61, 0xa3, 0x96, 0x17, 0x7a, 0x9c,
        0xb4, 0x10, 0xff, 0x61, 0xf2, 0x00, 0x15, 0xad,
    ];
    assert_eq!(abc_hash, expected, "SHA-256 of 'abc' should match known value");
}

/// Test Merkle tree pair hashing (internal node construction).
#[test]
fn test_merkle_pair_hash() {
    let left = [0u8; 32];
    let right = [1u8; 32];

    // Use SHA-256 explicitly for spec compatibility
    let result = Sha256Hash::hash_pair(&left, &right);
    assert_eq!(result.len(), 32, "Pair hash should be 32 bytes");

    // Verify it's not just concatenation
    assert_ne!(&result[..], &left[..], "Pair hash should not equal left");
    assert_ne!(&result[..], &right[..], "Pair hash should not equal right");

    // Verify determinism
    let result2 = Sha256Hash::hash_pair(&left, &right);
    assert_eq!(result, result2, "Pair hash should be deterministic");

    // Verify order matters
    let result_swapped = Sha256Hash::hash_pair(&right, &left);
    assert_ne!(result, result_swapped, "Pair hash should be order-sensitive");

    // Also test with the default hash_pair function
    let default_result = hash_pair(&left, &right);
    assert_eq!(default_result.len(), 32, "Default pair hash should be 32 bytes");
}

/// Test Merkle tree with known test vector.
/// This uses the test vector from the spec (SHA-256).
#[test]
fn test_merkle_tree_known_vector() {
    // Test vector from the Longfellow spec (requires SHA-256)
    let leaves: Vec<[u8; 32]> = vec![
        hex::decode("4bf5122f344554c53bde2ebb8cd2b7e3d1600ad631c385a5d7cce23c7785459a")
            .unwrap().try_into().unwrap(),
        hex::decode("dbc1b4c900ffe48d575b5da5c638040125f65db0fe3e24494b76ea986457d986")
            .unwrap().try_into().unwrap(),
        hex::decode("084fed08b978af4d7d196a7446a86b58009e636b611db16211b65a9aadff29c5")
            .unwrap().try_into().unwrap(),
        hex::decode("e52d9c508c502347344d8c07ad91cbd6068afc75ff6292f062a09ca381c89e71")
            .unwrap().try_into().unwrap(),
        hex::decode("e77b9a9ae9e30b0dbdb6f510a264ef9de781501d7b6b92ae89eb059c5ab743db")
            .unwrap().try_into().unwrap(),
    ];

    // Use SHA-256 explicitly to match spec test vector
    let mut tree = MerkleTreeGeneric::<Sha256Hash>::new(5);
    for (i, leaf) in leaves.iter().enumerate() {
        tree.set_leaf(i, *leaf);
    }

    let root = tree.build();

    // Expected root from spec test vector
    let expected_root = hex::decode("f22f4501ffd3bdffcecc9e4cd6828a4479aeedd6aa484eb7c1f808ccf71c6e76")
        .expect("valid hex");

    assert_eq!(
        &root[..], &expected_root[..],
        "Merkle root should match known test vector"
    );
}

/// Test polynomial evaluation at known points.
#[test]
fn test_polynomial_evaluation() {
    // Polynomial: p(x) = 1 + 2x + 3x^2
    let coeffs: Vec<Fp128> = vec![
        Fp128::from_u64(1),
        Fp128::from_u64(2),
        Fp128::from_u64(3),
    ];
    let poly = Polynomial::from_coeffs(coeffs);

    // p(0) = 1
    let p0 = poly.evaluate(Fp128::ZERO);
    assert_eq!(p0, Fp128::from_u64(1), "p(0) should equal 1");

    // p(1) = 1 + 2 + 3 = 6
    let p1 = poly.evaluate(Fp128::ONE);
    assert_eq!(p1, Fp128::from_u64(6), "p(1) should equal 6");

    // p(2) = 1 + 4 + 12 = 17
    let p2 = poly.evaluate(Fp128::from_u64(2));
    assert_eq!(p2, Fp128::from_u64(17), "p(2) should equal 17");
}

/// Test Lagrange interpolation roundtrip.
#[test]
fn test_lagrange_interpolation_roundtrip() {
    // Points: (0, 5), (1, 7), (2, 13)
    // These fit the polynomial p(x) = 5 + x + x^2
    let x_vals: Vec<Fp128> = vec![
        Fp128::from_u64(0),
        Fp128::from_u64(1),
        Fp128::from_u64(2),
    ];
    let y_vals: Vec<Fp128> = vec![
        Fp128::from_u64(5),
        Fp128::from_u64(7),
        Fp128::from_u64(13),
    ];

    let poly = Polynomial::interpolate(&x_vals, &y_vals);

    // Verify: evaluating at original points should give original values
    for (x, y) in x_vals.iter().zip(y_vals.iter()) {
        let eval = poly.evaluate(*x);
        assert_eq!(eval, *y, "Interpolation should recover original points");
    }
}

/// Test polynomial extension (Reed-Solomon encoding).
#[test]
fn test_polynomial_extension() {
    // Simple test: extend [1, 2, 3] from block=3 to ncol=6
    let coeffs: Vec<Fp128> = vec![
        Fp128::from_u64(1),
        Fp128::from_u64(2),
        Fp128::from_u64(3),
    ];

    let extended = extend(&coeffs, 3, 6);
    assert_eq!(extended.len(), 6, "Extended vector should have length ncol");

    // The first 3 elements should be evaluations at 0, 1, 2
    // For coeffs [1, 2, 3]: p(x) = 1 + 2x + 3x^2
    // p(0) = 1, p(1) = 6, p(2) = 17
    assert_eq!(extended[0], Fp128::from_u64(1), "p(0) = 1");
    assert_eq!(extended[1], Fp128::from_u64(6), "p(1) = 6");
    assert_eq!(extended[2], Fp128::from_u64(17), "p(2) = 17");
}

/// Test transcript determinism.
#[test]
fn test_transcript_determinism() {
    let mut t1 = Transcript::new(b"test-domain");
    let mut t2 = Transcript::new(b"test-domain");

    // Write same data to both
    t1.write_bytes(b"hello");
    t2.write_bytes(b"hello");

    // Generate challenges
    let c1: Vec<Fp128> = t1.generate_challenges(3);
    let c2: Vec<Fp128> = t2.generate_challenges(3);

    assert_eq!(c1, c2, "Transcripts with same input should produce same challenges");
}

/// Test that different domain separators produce different challenges.
#[test]
fn test_transcript_domain_separation() {
    let mut t1 = Transcript::new(b"domain-a");
    let mut t2 = Transcript::new(b"domain-b");

    let c1: Vec<Fp128> = t1.generate_challenges(1);
    let c2: Vec<Fp128> = t2.generate_challenges(1);

    assert_ne!(c1, c2, "Different domains should produce different challenges");
}

/// Test transcript index generation without replacement.
#[test]
fn test_transcript_indices_without_replacement() {
    let mut transcript = Transcript::new(b"test");
    transcript.write_bytes(b"seed");

    let indices = transcript.generate_indices_without_replacement(100, 10);

    // Should have exactly 10 indices
    assert_eq!(indices.len(), 10, "Should generate 10 indices");

    // All indices should be unique
    let mut sorted = indices.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(sorted.len(), 10, "All indices should be unique");

    // All indices should be in range [0, 100)
    for &idx in &indices {
        assert!(idx < 100, "Index should be in range");
    }
}

/// Test Ligero parameters match specification.
#[test]
fn test_ligero_params_spec() {
    // For 128-bit security, test that parameters are reasonable
    let params = LigeroParams::new(100, 10, 128);

    // NREQ should be around 6 for 128-bit security
    assert!(params.nreq >= 4, "NREQ should be at least 4 for security");
    assert!(params.nreq <= 12, "NREQ should not be excessive");

    // RATE should be 4 (standard value)
    assert_eq!(params.rate, 4, "Rate should be 4");

    // WR >= NREQ
    assert!(params.wr >= params.nreq, "WR should be >= NREQ");

    // BLOCK = NREQ + WR
    assert!(params.block >= params.nreq + params.wr, "BLOCK should be >= NREQ + WR");

    // DBLOCK = 2 * BLOCK - 1
    assert_eq!(params.dblock, 2 * params.block - 1, "DBLOCK should be 2*BLOCK - 1");

    // NCOL = BLOCK * RATE
    assert_eq!(params.ncol, params.block * params.rate, "NCOL should be BLOCK * RATE");
}

/// Test proof serialization size is reasonable.
#[test]
fn test_ligero_proof_size() {
    use longfellow_zk::merkle::MerkleProof;

    // Create a minimal proof structure
    let ldt: Vec<Fp128> = vec![Fp128::from_u64(1), Fp128::from_u64(2)];
    let dot: Vec<Fp128> = vec![Fp128::from_u64(3), Fp128::from_u64(4)];
    let qpr: Vec<Fp128> = vec![Fp128::from_u64(5), Fp128::from_u64(6)];
    let columns: Vec<Vec<Fp128>> = vec![
        vec![Fp128::from_u64(7), Fp128::from_u64(8)],
        vec![Fp128::from_u64(9), Fp128::from_u64(10)],
    ];

    // Create empty merkle proof from bytes
    let merkle_bytes: Vec<u8> = vec![
        5, 0, 0, 0, 0, 0, 0, 0,  // leaf_count = 5
        0, 0, 0, 0, 0, 0, 0, 0,  // node_count = 0
    ];
    let merkle_proof = MerkleProof::from_bytes(&merkle_bytes).expect("valid proof bytes");

    let proof = LigeroProof::new(ldt, dot, qpr, columns, merkle_proof);

    // Serialize
    let bytes = proof.to_bytes();

    // Verify structure: should have length prefixes
    assert!(!bytes.is_empty(), "Serialized proof should not be empty");

    // Verify size method matches
    assert_eq!(proof.size(), bytes.len(), "Size method should match actual size");
}

/// Test that field element distribution in challenges appears random.
#[test]
fn test_challenge_distribution() {
    let mut transcript = Transcript::new(b"distribution-test");

    // Generate many challenges
    let challenges: Vec<Fp128> = transcript.generate_challenges(100);

    // Basic sanity checks
    assert_eq!(challenges.len(), 100, "Should generate 100 challenges");

    // No two consecutive challenges should be equal (with high probability)
    let mut duplicates = 0;
    for i in 1..challenges.len() {
        if challenges[i] == challenges[i - 1] {
            duplicates += 1;
        }
    }
    assert!(duplicates < 5, "Should have very few consecutive duplicates");
}

/// Test large field elements near the modulus boundary.
#[test]
fn test_field_boundary_values() {
    // Test values near 2^64
    let two_64 = Fp128::from_raw([0, 1]);
    let two_64_minus_1 = Fp128::from_raw([u64::MAX, 0]);

    let sum = two_64 + two_64_minus_1;
    let expected = Fp128::from_raw([u64::MAX, 1]);
    assert_eq!(sum, expected, "2^64 + (2^64 - 1) should equal 2^65 - 1");

    // Test values near 2^108: 2^108 = 0x0000100000000000_0000000000000000
    let two_108 = Fp128::from_raw([0, 0x0000100000000000]);

    // 2^108 * 2 = 2^109
    let two_109 = two_108 + two_108;
    let expected_109 = Fp128::from_raw([0, 0x0000200000000000]);
    assert_eq!(two_109, expected_109, "2^108 + 2^108 should equal 2^109");
}

/// Test multiplication by powers of 2.
#[test]
fn test_field_powers_of_two() {
    let two = Fp128::from_u64(2);
    let mut power = Fp128::ONE;

    // Compute 2^10 = 1024
    for _ in 0..10 {
        power = power * two;
    }
    assert_eq!(power, Fp128::from_u64(1024), "2^10 should equal 1024");

    // Compute 2^20 = 1048576
    for _ in 0..10 {
        power = power * two;
    }
    assert_eq!(power, Fp128::from_u64(1048576), "2^20 should equal 1048576");
}

/// Test polynomial interpolation with more points.
#[test]
fn test_interpolation_larger() {
    // Test with 10 points
    let n = 10;
    let x_vals: Vec<Fp128> = (0..n).map(|i| Fp128::from_u64(i as u64)).collect();
    let y_vals: Vec<Fp128> = (0..n).map(|i| Fp128::from_u64((i * i + 1) as u64)).collect();

    let poly = Polynomial::interpolate(&x_vals, &y_vals);

    // Verify at original points
    for (x, y) in x_vals.iter().zip(y_vals.iter()) {
        let eval = poly.evaluate(*x);
        assert_eq!(eval, *y, "Interpolation should recover original points");
    }
}

/// Test end-to-end ZK proof generation and verification.
#[test]
fn test_zk_proof_roundtrip() {
    use longfellow_zk::circuit::{CircuitBuilder, LayerBuilder};
    use longfellow_zk::zk::{verify_zk, ZkProver};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    // Create a simple circuit: a * b = c
    let mut builder = CircuitBuilder::<Fp128>::new(0, 3);
    let mut layer = LayerBuilder::new(0, 2);
    layer.add_mul(0, 0, 1, Fp128::ONE);
    builder.add_layer(layer);
    let circuit = builder.build();

    // Create witness: 3 * 4 = 12
    let witness = vec![
        Fp128::from_u64(3),
        Fp128::from_u64(4),
        Fp128::from_u64(12),
    ];

    // Create prover and generate proof
    let prover = ZkProver::new(circuit.clone(), vec![], witness);
    let mut rng = ChaCha20Rng::seed_from_u64(12345);
    let proof = prover.prove(&mut rng);

    // Verify proof
    let verified = verify_zk(&circuit, &[], &proof);
    assert!(verified, "Valid proof should verify");
}

// ==========================================================================
// IETF Spec Field Compatibility Tests
// ==========================================================================

/// Test GF(2^128) irreducible polynomial matches the IETF spec.
/// Spec: GF(2)[x] / (x^128 + x^7 + x^2 + x + 1)
#[test]
fn test_gf2_128_irreducible_polynomial() {
    // The defining relation: x^128 = x^7 + x^2 + x + 1 = 0x87
    let x = GF2_128::X;
    let x128 = x.pow(128);
    assert_eq!(x128, GF2_128::from_u64(0x87),
        "x^128 should reduce to x^7 + x^2 + x + 1 per IETF spec");

    // Verify x is a generator of the multiplicative group
    // (order of group is 2^128 - 1, so x^(2^128 - 1) = 1)
    // We don't test the full order, but x should not be 0 or 1
    assert_ne!(x, GF2_128::ZERO);
    assert_ne!(x, GF2_128::ONE);
}

/// Test GF(2^128) P2 constant for sumcheck matches spec.
/// Spec: For GF(2^128), P2 = X (the polynomial variable).
#[test]
fn test_gf2_128_p2_constant() {
    // P2 = X per spec (NOT 2, since 1+1=0 in characteristic 2)
    assert_eq!(<GF2_128 as Field>::P2, GF2_128::X);

    // Verify P2 is distinct from P0=0 and P1=1
    assert_ne!(<GF2_128 as Field>::P2, GF2_128::ZERO);
    assert_ne!(<GF2_128 as Field>::P2, GF2_128::ONE);
}

/// Test GF(2^128) characteristic-2 properties.
#[test]
fn test_gf2_128_characteristic_2() {
    // a + a = 0 for all a
    let a = GF2_128::from_raw([0xdeadbeef, 0xcafebabe]);
    assert_eq!(a + a, GF2_128::ZERO);

    // -a = a
    assert_eq!(-a, a);

    // Frobenius: (a+b)^2 = a^2 + b^2
    let b = GF2_128::from_raw([0x12345678, 0x9abcdef0]);
    assert_eq!((a + b).square(), a.square() + b.square());
}

/// Test GF(2^128) serialization: 16 bytes, little-endian polynomial.
/// Spec: "obtain b = fs.bytes(16) and interpret the 128 bits of b as a
/// little-endian polynomial"
#[test]
fn test_gf2_128_serialization_spec() {
    // X = polynomial variable = bit 1 set = byte[0] = 0x02
    let x = GF2_128::X;
    let bytes = x.to_bytes();
    assert_eq!(bytes.len(), 16);
    assert_eq!(bytes[0], 0x02);
    for b in &bytes[1..] {
        assert_eq!(*b, 0);
    }

    // x^7 + x^2 + x + 1 = 0x87 = byte[0] = 0x87
    let poly = GF2_128::from_u64(0x87);
    let bytes = poly.to_bytes();
    assert_eq!(bytes[0], 0x87);
    for b in &bytes[1..] {
        assert_eq!(*b, 0);
    }
}

/// Test Fp256 modulus matches NIST P-256 base field.
/// Spec: p = 2^256 - 2^224 + 2^192 + 2^96 - 1
/// = 115792089210356248762697446949407573530086143415290314195533631308867097853951
/// = 0xffffffff00000001_0000000000000000_00000000ffffffff_ffffffffffffffff
#[test]
fn test_fp256_modulus_spec() {
    // Compute p - 1 and verify (p - 1) + 1 = 0 mod p
    let p_minus_1 = Fp256::ZERO - Fp256::ONE;
    assert_eq!(p_minus_1 + Fp256::ONE, Fp256::ZERO);

    // Verify the specific byte pattern of p-1 (little-endian serialization)
    // p - 1 = 0xffffffff00000001_0000000000000000_00000000ffffffff_fffffffffffffffe
    let bytes = p_minus_1.to_bytes();
    // Byte 0 (LSB): 0xfe (since p ends in ...ffff and p-1 ends in ...fffe)
    assert_eq!(bytes[0], 0xfe, "p-1 LSB should be 0xfe");
    // Bytes 1..8 (rest of limb[0]): all 0xff
    assert_eq!(&bytes[1..8], &[0xff; 7], "limb[0] high bytes");
    // Bytes 8..12 (low half of limb[1] = 0xffffffff): all 0xff
    assert_eq!(&bytes[8..12], &[0xff; 4], "limb[1] low half");
    // Bytes 12..16 (high half of limb[1] = 0x00000000): all 0x00
    assert_eq!(&bytes[12..16], &[0x00; 4], "limb[1] high half");
    // Bytes 16..24 (limb[2] = 0): all 0x00
    assert_eq!(&bytes[16..24], &[0x00; 8], "limb[2]");
    // Bytes 24..28 (low half of limb[3] = 0x00000001): 01 00 00 00
    assert_eq!(&bytes[24..28], &[0x01, 0x00, 0x00, 0x00], "limb[3] low half");
    // Bytes 28..32 (high half of limb[3] = 0xffffffff): all 0xff
    assert_eq!(&bytes[28..32], &[0xff; 4], "limb[3] high half");
}

/// Test Fp256 serialization: 32 bytes, little-endian.
/// Spec: "a 256-bit element of the finite field Fp256 is serialized into
/// 32-bytes starting with the least-significant byte"
#[test]
fn test_fp256_serialization_spec() {
    let val = Fp256::from_u64(0x1234);
    let bytes = val.to_bytes();
    assert_eq!(bytes.len(), 32);
    // Little-endian: low bytes first
    assert_eq!(bytes[0], 0x34);
    assert_eq!(bytes[1], 0x12);
    for b in &bytes[2..] {
        assert_eq!(*b, 0);
    }

    // Round-trip
    let recovered = Fp256::from_bytes(&bytes);
    assert_eq!(recovered, val);
}

/// Test FieldID values match the IETF spec assignments.
/// Spec assigns: p256=0x01, p384=0x02, p521=0x03, GF(2^128)=0x04,
/// GF(2^16)=0x05, 2^128-2^108+1=0x06, 2^64-59=0x07, 2^64-2^32+1=0x08,
/// F_{2^64-59}^2=0x09, secp256=0x0a
#[test]
fn test_field_id_spec_values() {
    assert_eq!(FieldId::P256 as u8, 0x01);
    assert_eq!(FieldId::P384 as u8, 0x02);
    assert_eq!(FieldId::P521 as u8, 0x03);
    assert_eq!(FieldId::Gf2_128 as u8, 0x04);
    assert_eq!(FieldId::Gf2_16 as u8, 0x05);
    assert_eq!(FieldId::Fp128 as u8, 0x06);
    assert_eq!(FieldId::Fp64m59 as u8, 0x07);
    assert_eq!(FieldId::Goldilocks as u8, 0x08);
    assert_eq!(FieldId::Fp64m59Sq as u8, 0x09);
    assert_eq!(FieldId::Secp256 as u8, 0x0a);

    // Verify implemented fields have correct IDs
    assert_eq!(Fp256::FIELD_ID, FieldId::P256);
    assert_eq!(GF2_128::FIELD_ID, FieldId::Gf2_128);
    assert_eq!(Fp128::FIELD_ID, FieldId::Fp128);
}

/// Test GF(2^128) x^(-1) = x^127 + x^6 + x + 1 (from reference implementation).
#[test]
fn test_gf2_128_x_inverse_known_value() {
    // From google/longfellow-zk: x^{-1} = x^127 + x^6 + x + 1
    // = 0x8000000000000000_0000000000000043 (hi = 2^63, lo = 0x43 = 0b1000011)
    assert_eq!(GF2_128::X * GF2_128::X_INV, GF2_128::ONE);

    // Also verify via general inversion
    let x_inv_computed = GF2_128::X.invert().unwrap();
    assert_eq!(x_inv_computed, GF2_128::X_INV);
}

/// Test that ZK proofs work end-to-end over all spec-required fields.
#[test]
fn test_zk_all_fields() {
    use longfellow_zk::circuit::{CircuitBuilder, LayerBuilder};
    use longfellow_zk::zk::{verify_zk, ZkProver};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    fn test_field<F: Field>() -> bool {
        let mut builder = CircuitBuilder::<F>::new(0, 3);
        let mut layer = LayerBuilder::new(0, 2);
        layer.add_mul(0, 0, 1, F::ONE);
        builder.add_layer(layer);
        let circuit = builder.build();

        let a = F::from_u64(5);
        let b = F::from_u64(9);
        let witness = vec![a, b, a * b];

        let prover = ZkProver::new(circuit.clone(), vec![], witness);
        let mut rng = ChaCha20Rng::seed_from_u64(777);
        let proof = prover.prove(&mut rng);

        verify_zk(&circuit, &[], &proof)
    }

    assert!(test_field::<Fp128>(), "Fp128 ZK proof should verify");
    assert!(test_field::<GF2_128>(), "GF(2^128) ZK proof should verify");
    assert!(test_field::<Fp256>(), "Fp256 ZK proof should verify");
}
