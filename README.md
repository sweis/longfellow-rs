# Longfellow ZK

A Rust implementation of the Longfellow zero-knowledge proof system based on the IETF draft specification.

⚠️ This is experimental and written by Claude. It has had almost no review. Do not use for real applications yet.

## Overview

Longfellow is a succinct non-interactive zero-knowledge (SNARK) proof system that combines:
- **Ligero commitment scheme** - A polynomial commitment based on Reed-Solomon codes and Merkle trees
- **Sumcheck protocol** - For verifying polynomial identities over layered arithmetic circuits

The scheme provides 128-bit security with efficient prover and verifier algorithms.

## References

### Specification
- [IETF Draft: draft-google-cfrg-libzk-01](https://datatracker.ietf.org/doc/draft-google-cfrg-libzk/)

### Reference Implementations
- [google/longfellow-zk](https://github.com/google/longfellow-zk) - C++ reference implementation
- [abetterinternet/zk-cred-longfellow](https://github.com/abetterinternet/zk-cred-longfellow) - Rust implementation

### Papers
- [Ligero: Lightweight Sublinear Arguments Without a Trusted Setup](https://eprint.iacr.org/2022/1608) - The underlying commitment scheme
- [Proofs, Arguments, and Zero-Knowledge](https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf) - Background on sumcheck protocol

## Building

### Prerequisites
- Rust 1.70 or later
- Cargo

### Build
```bash
cargo build --release
```

### Run Tests
```bash
cargo test
```

### Run Benchmarks
```bash
cargo bench
```

## Usage

Add to your `Cargo.toml`:
```toml
[dependencies]
longfellow-zk = { path = "." }
```

### Examples

| Example | Description | Command |
|---------|-------------|---------|
| `simple_proof` | Basic multiplication circuit | `cargo run --example simple_proof` |
| `advanced_circuits` | Age verification, inner product, etc. | `cargo run --example advanced_circuits` |
| `multi_field` | Proofs over Fp128, GF(2^128), and Fp256 | `cargo run --example multi_field` |

## Architecture

```
src/
├── lib.rs          # Public API and module exports
├── field/          # Finite field arithmetic
│   ├── mod.rs      # Field trait + FieldId + batch_invert
│   ├── fp128.rs    # Prime field 2^128 - 2^108 + 1 (FieldID 0x06)
│   ├── fp256.rs    # NIST P-256 base field (FieldID 0x01)
│   └── gf2_128.rs  # Binary extension GF(2^128) (FieldID 0x04)
├── polynomial.rs   # Polynomial operations and Reed-Solomon extension
├── merkle.rs       # Merkle tree for commitments
├── transcript.rs   # Fiat-Shamir transcript (SHA-256 based)
├── circuit.rs      # Circuit representation
├── zk.rs           # High-level ZK prover and verifier
├── ligero/         # Ligero commitment scheme
│   ├── mod.rs
│   ├── params.rs   # Ligero parameters
│   ├── commitment.rs
│   ├── prover.rs   # Ligero prover
│   ├── verifier.rs # Ligero verifier
│   └── proof.rs    # Proof structure
└── sumcheck/       # Sumcheck protocol
    ├── mod.rs
    ├── eq.rs       # Equality polynomial operations
    ├── layer.rs    # Layer-by-layer sumcheck prover
    ├── proof.rs    # Full circuit sumcheck prover
    └── verifier.rs # Sumcheck verifier
```

## Supported Fields

All fields required by the IETF draft-google-cfrg-libzk are implemented:

| Field     | FieldID | Characteristic | Definition                              | Use Case              |
|-----------|---------|----------------|-----------------------------------------|-----------------------|
| `Fp256`   | 0x01    | large prime    | p = 2^256 - 2^224 + 2^192 + 2^96 - 1    | ECDSA-P256 signatures |
| `GF2_128` | 0x04    | 2              | GF(2)[x] / (x^128 + x^7 + x^2 + x + 1)  | SHA-256 circuits      |
| `Fp128`   | 0x06    | large prime    | p = 2^128 - 2^108 + 1                   | General arithmetic    |

The entire proof system (Ligero, sumcheck, circuits) is generic over any
type implementing the `Field` trait.

### Binary Field Notes

GF(2^128) has characteristic 2, which means:
- Addition is XOR: `a + a = 0`
- Negation is identity: `-a = a`
- The sumcheck P2 constant is the polynomial variable X (not 2, since 1+1=0)
- Multiplication is carryless (polynomial multiplication mod Q(x))

## Security

- **Security level**: 128 bits
- **Hash function**: SHA-256 for Fiat-Shamir transform and Merkle tree
- **Soundness**: Based on Reed-Solomon proximity testing

## Testing

Run the full test suite:
```bash
cargo test
```

Run specific test modules:
```bash
cargo test field::             # All field arithmetic tests (Fp128, Fp256, GF2_128)
cargo test polynomial::tests   # Polynomial operations
cargo test merkle::tests       # Merkle tree tests
cargo test ligero::            # Ligero commitment tests
cargo test sumcheck::          # Sumcheck protocol tests (prover + verifier)
cargo test zk::tests           # End-to-end ZK tests over all fields
```

## Verification

The `verify_zk` function performs complete end-to-end verification:

1. **Sumcheck verification**: replays the Fiat-Shamir transcript round-by-round,
   checking that each polynomial evaluation reduces the claim correctly.
2. **Ligero verification**:
   - Merkle proof check (committed columns match the root)
   - Low-degree test (each row is a valid polynomial)
   - Dot-product test (linear constraints on witness)
   - Quadratic test (x·y = z constraints)

Tampered proofs (modified commitment, LDT response, or layer count) are rejected.
