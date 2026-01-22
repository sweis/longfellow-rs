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

### Basic Example

See `examples/simple_proof.rs` for a complete example demonstrating:
- Creating a circuit with quadratic constraints
- Generating a zero-knowledge proof
- Verifying the proof

Run the example:
```bash
cargo run --example simple_proof
```

## Architecture

```
src/
├── lib.rs          # Public API and module exports
├── field.rs        # Fp128 finite field arithmetic (p = 2^128 - 2^108 + 1)
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
    ├── layer.rs    # Layer-by-layer sumcheck
    └── proof.rs    # Sumcheck prover
```

## Field Arithmetic

The implementation uses the prime field Fp128 with:
- Prime: p = 2^128 - 2^108 + 1
- This prime allows efficient reduction using the identity: 2^128 ≡ 2^108 - 1 (mod p)

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
cargo test field::tests        # Field arithmetic tests
cargo test polynomial::tests   # Polynomial operations
cargo test merkle::tests       # Merkle tree tests
cargo test ligero::tests       # Ligero commitment tests
cargo test sumcheck::tests     # Sumcheck protocol tests
cargo test zk::tests           # End-to-end ZK tests
```

## Contributing

Contributions are welcome! Please ensure all tests pass before submitting a pull request.
