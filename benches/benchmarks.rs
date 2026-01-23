//! Benchmarks for the Longfellow ZK library.
//!
//! Run with: cargo bench
//! Run specific benchmark: cargo bench -- "Field"

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use longfellow_zk::{
    batch_invert,
    circuit::{CircuitBuilder, LayerBuilder},
    fft::{fft, ifft, polynomial_multiply, FftDomain},
    field::Fp128,
    hash::{HashFunction, Sha256Hash},
    merkle::{MerkleTree, MerkleTreeGeneric},
    polynomial::{extend, extend_evaluations, Polynomial},
    transcript::Transcript,
    zk::{verify_zk, ZkProver},
    ligero::LigeroProver,
};
#[cfg(feature = "blake3_hash")]
use longfellow_zk::hash::Blake3Hash;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

fn bench_field_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("Field Operations");

    let a = Fp128::from_u64(12345678901234567890);
    let b = Fp128::from_u64(9876543210987654321);

    group.bench_function("add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
    });

    group.bench_function("sub", |bencher| {
        bencher.iter(|| black_box(a) - black_box(b))
    });

    group.bench_function("mul", |bencher| {
        bencher.iter(|| black_box(a) * black_box(b))
    });

    group.bench_function("square", |bencher| {
        bencher.iter(|| black_box(a).square())
    });

    group.bench_function("invert", |bencher| {
        bencher.iter(|| black_box(a).invert())
    });

    group.bench_function("negate", |bencher| {
        bencher.iter(|| -black_box(a))
    });

    // Benchmark batch inversion
    for size in [10, 50, 100, 500].iter() {
        let elements: Vec<Fp128> = (1..=*size)
            .map(|i| Fp128::from_u64(i as u64))
            .collect();

        group.bench_with_input(BenchmarkId::new("batch_invert", size), size, |b, _| {
            b.iter(|| black_box(batch_invert(&elements)))
        });

        // Compare with sequential inversion
        group.bench_with_input(BenchmarkId::new("sequential_invert", size), size, |b, _| {
            b.iter(|| {
                let result: Vec<_> = elements.iter().map(|e| e.invert().unwrap()).collect();
                black_box(result)
            })
        });
    }

    group.finish();
}

fn bench_polynomial(c: &mut Criterion) {
    let mut group = c.benchmark_group("Polynomial");

    // Test polynomial extension at different sizes
    for size in [16, 32, 64, 128].iter() {
        let coeffs: Vec<Fp128> = (0..*size)
            .map(|i| Fp128::from_u64(i as u64 + 1))
            .collect();

        group.throughput(Throughput::Elements(*size as u64));

        group.bench_with_input(BenchmarkId::new("extend", size), size, |b, &size| {
            b.iter(|| black_box(extend(&coeffs, size, size * 4)))
        });

        group.bench_with_input(BenchmarkId::new("extend_evaluations", size), size, |b, &size| {
            let evals = extend(&coeffs, size, size);
            b.iter(|| black_box(extend_evaluations(&evals, size, size * 4)))
        });
    }

    // Test polynomial interpolation at different sizes
    for size in [8, 16, 32].iter() {
        let x_vals: Vec<Fp128> = (0..*size)
            .map(|i| Fp128::from_u64(i as u64))
            .collect();
        let y_vals: Vec<Fp128> = (0..*size)
            .map(|i| Fp128::from_u64((i * i + 1) as u64))
            .collect();

        group.bench_with_input(BenchmarkId::new("interpolate", size), size, |b, _| {
            b.iter(|| black_box(Polynomial::interpolate(&x_vals, &y_vals)))
        });
    }

    // Test polynomial evaluation
    let coeffs: Vec<Fp128> = (0..32)
        .map(|i| Fp128::from_u64(i as u64 + 1))
        .collect();
    let poly = Polynomial::from_coeffs(coeffs);
    let x = Fp128::from_u64(42);

    group.bench_function("evaluate_32", |bencher| {
        bencher.iter(|| black_box(poly.evaluate(black_box(x))))
    });

    group.finish();
}

fn bench_merkle_tree(c: &mut Criterion) {
    let mut group = c.benchmark_group("Merkle Tree");

    for size in [16, 64, 256, 1024].iter() {
        group.throughput(Throughput::Elements(*size as u64));

        group.bench_with_input(BenchmarkId::new("build", size), size, |b, &size| {
            let mut tree = MerkleTree::new(size);
            for i in 0..size {
                let mut leaf = [0u8; 32];
                leaf[0..8].copy_from_slice(&(i as u64).to_le_bytes());
                tree.set_leaf(i, leaf);
            }
            b.iter(|| {
                let mut t = tree.clone();
                black_box(t.build())
            });
        });

        group.bench_with_input(BenchmarkId::new("proof", size), size, |b, &size| {
            let mut tree = MerkleTree::new(size);
            for i in 0..size {
                let mut leaf = [0u8; 32];
                leaf[0..8].copy_from_slice(&(i as u64).to_le_bytes());
                tree.set_leaf(i, leaf);
            }
            tree.build();

            // Request proof for first 6 indices
            let indices: Vec<usize> = (0..6.min(size)).collect();
            b.iter(|| {
                black_box(tree.compressed_proof(&indices))
            });
        });
    }

    group.finish();
}

fn bench_transcript(c: &mut Criterion) {
    let mut group = c.benchmark_group("Transcript");

    group.bench_function("generate_challenge", |bencher| {
        let mut transcript = Transcript::new(b"benchmark");
        transcript.write_bytes(&[0u8; 32]);
        bencher.iter(|| {
            let c: Fp128 = transcript.generate_challenge();
            black_box(c)
        })
    });

    group.bench_function("generate_challenges_10", |bencher| {
        bencher.iter(|| {
            let mut transcript = Transcript::new(b"benchmark");
            transcript.write_bytes(&[0u8; 32]);
            let c: Vec<Fp128> = transcript.generate_challenges(10);
            black_box(c)
        })
    });

    group.bench_function("write_field_element", |bencher| {
        let elem = Fp128::from_u64(42);
        bencher.iter(|| {
            let mut transcript = Transcript::new(b"benchmark");
            transcript.write_field_element(&elem);
            black_box(transcript)
        })
    });

    group.finish();
}

fn bench_ligero(c: &mut Criterion) {
    let mut group = c.benchmark_group("Ligero");
    group.sample_size(20); // Reduce sample size for slower benchmarks

    // Test with different witness sizes
    for witness_size in [10, 50, 100].iter() {
        let witness: Vec<Fp128> = (0..*witness_size)
            .map(|i| Fp128::from_u64(i as u64 + 1))
            .collect();

        group.throughput(Throughput::Elements(*witness_size as u64));

        group.bench_with_input(BenchmarkId::new("commit", witness_size), witness_size, |b, _| {
            let mut rng = ChaCha20Rng::seed_from_u64(42);
            b.iter(|| {
                let mut prover = LigeroProver::<Fp128>::new(witness.clone(), vec![], 128);
                black_box(prover.commit(&mut rng))
            });
        });
    }

    group.finish();
}

fn bench_simple_proof(c: &mut Criterion) {
    let mut group = c.benchmark_group("Simple Proof");
    group.sample_size(20);

    // Create a simple multiplication circuit
    let mut builder = CircuitBuilder::<Fp128>::new(0, 3);
    let mut layer = LayerBuilder::new(0, 2);
    layer.add_mul(0, 0, 1, Fp128::ONE);
    builder.add_layer(layer);
    let circuit = builder.build();

    let witness = vec![
        Fp128::from_u64(3),
        Fp128::from_u64(4),
        Fp128::from_u64(12),
    ];

    group.bench_function("prove", |bencher| {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        bencher.iter(|| {
            let prover = ZkProver::new(circuit.clone(), vec![], witness.clone());
            black_box(prover.prove(&mut rng))
        })
    });

    group.bench_function("verify", |bencher| {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let prover = ZkProver::new(circuit.clone(), vec![], witness.clone());
        let proof = prover.prove(&mut rng);

        bencher.iter(|| {
            black_box(verify_zk(&circuit, &[], &proof))
        })
    });

    group.finish();
}

fn bench_larger_circuit(c: &mut Criterion) {
    let mut group = c.benchmark_group("Larger Circuit");
    group.sample_size(10); // Reduce for slower benchmarks

    // Create a circuit with multiple multiplication gates
    // Circuit: x0 * x1 = x2, x2 * x3 = x4, x4 * x5 = x6, etc.
    for num_muls in [2, 5, 10].iter() {
        let witness_size = num_muls * 2 + 1;

        let mut builder = CircuitBuilder::<Fp128>::new(0, witness_size);
        for i in 0..*num_muls {
            let mut layer = LayerBuilder::new(0, 2);
            layer.add_mul(0, i * 2, i * 2 + 1, Fp128::ONE);
            builder.add_layer(layer);
        }
        let circuit = builder.build();

        // Create witness: each pair multiplies to the next value
        let mut witness = Vec::with_capacity(witness_size);
        witness.push(Fp128::from_u64(2));
        for i in 0..*num_muls {
            witness.push(Fp128::from_u64((i + 2) as u64));
            // Product of previous two
            let product = witness[i * 2] * witness[i * 2 + 1];
            if i < *num_muls - 1 {
                witness.push(product);
            }
        }
        // Ensure witness is correct length
        while witness.len() < witness_size {
            witness.push(Fp128::from_u64(1));
        }

        group.bench_with_input(BenchmarkId::new("prove", num_muls), num_muls, |b, _| {
            let mut rng = ChaCha20Rng::seed_from_u64(42);
            b.iter(|| {
                let prover = ZkProver::new(circuit.clone(), vec![], witness.clone());
                black_box(prover.prove(&mut rng))
            });
        });

        group.bench_with_input(BenchmarkId::new("verify", num_muls), num_muls, |b, _| {
            let mut rng = ChaCha20Rng::seed_from_u64(42);
            let prover = ZkProver::new(circuit.clone(), vec![], witness.clone());
            let proof = prover.prove(&mut rng);
            b.iter(|| {
                black_box(verify_zk(&circuit, &[], &proof))
            });
        });
    }

    group.finish();
}

fn bench_hash_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("Hash Comparison");

    // Test data of various sizes
    let data_small = [0u8; 64];
    let data_medium = [0u8; 256];
    let data_large = [0u8; 1024];

    // SHA-256 benchmarks
    group.bench_function("SHA256/hash_64B", |b| {
        b.iter(|| black_box(Sha256Hash::hash(&data_small)))
    });

    group.bench_function("SHA256/hash_256B", |b| {
        b.iter(|| black_box(Sha256Hash::hash(&data_medium)))
    });

    group.bench_function("SHA256/hash_1KB", |b| {
        b.iter(|| black_box(Sha256Hash::hash(&data_large)))
    });

    // SHA-256 hash_pair
    let left = [1u8; 32];
    let right = [2u8; 32];
    group.bench_function("SHA256/hash_pair", |b| {
        b.iter(|| black_box(Sha256Hash::hash_pair(&left, &right)))
    });

    // BLAKE3 benchmarks (when feature enabled)
    #[cfg(feature = "blake3_hash")]
    {
        group.bench_function("BLAKE3/hash_64B", |b| {
            b.iter(|| black_box(Blake3Hash::hash(&data_small)))
        });

        group.bench_function("BLAKE3/hash_256B", |b| {
            b.iter(|| black_box(Blake3Hash::hash(&data_medium)))
        });

        group.bench_function("BLAKE3/hash_1KB", |b| {
            b.iter(|| black_box(Blake3Hash::hash(&data_large)))
        });

        group.bench_function("BLAKE3/hash_pair", |b| {
            b.iter(|| black_box(Blake3Hash::hash_pair(&left, &right)))
        });
    }

    group.finish();
}

fn bench_fft(c: &mut Criterion) {
    let mut group = c.benchmark_group("FFT");

    // Test FFT at different sizes
    for log_size in [6, 8, 10, 12].iter() {
        let size = 1 << log_size;
        let domain = FftDomain::new(size);

        // Random coefficients
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let coeffs: Vec<Fp128> = (0..size).map(|_| Fp128::random(&mut rng)).collect();

        group.throughput(Throughput::Elements(size as u64));

        group.bench_with_input(BenchmarkId::new("fft", size), &size, |b, _| {
            b.iter(|| black_box(fft(&coeffs, &domain)))
        });

        group.bench_with_input(BenchmarkId::new("ifft", size), &size, |b, _| {
            let evals = fft(&coeffs, &domain);
            b.iter(|| black_box(ifft(&evals, &domain)))
        });
    }

    // Benchmark polynomial multiplication via FFT vs naive
    for size in [32, 64, 128].iter() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let a: Vec<Fp128> = (0..*size).map(|_| Fp128::random(&mut rng)).collect();
        let b: Vec<Fp128> = (0..*size).map(|_| Fp128::random(&mut rng)).collect();

        group.bench_with_input(BenchmarkId::new("poly_mul_fft", size), size, |bench, _| {
            bench.iter(|| black_box(polynomial_multiply(&a, &b)))
        });
    }

    group.finish();
}

fn bench_merkle_hash_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("Merkle Hash Comparison");
    group.sample_size(20);

    for size in [64, 256, 1024].iter() {
        // SHA-256 Merkle tree
        group.bench_with_input(BenchmarkId::new("SHA256/build", size), size, |b, &size| {
            let mut tree = MerkleTreeGeneric::<Sha256Hash>::new(size);
            for i in 0..size {
                let mut leaf = [0u8; 32];
                leaf[0..8].copy_from_slice(&(i as u64).to_le_bytes());
                tree.set_leaf(i, leaf);
            }
            b.iter(|| {
                let mut t = tree.clone();
                black_box(t.build())
            });
        });

        // BLAKE3 Merkle tree (when feature enabled)
        #[cfg(feature = "blake3_hash")]
        group.bench_with_input(BenchmarkId::new("BLAKE3/build", size), size, |b, &size| {
            let mut tree = MerkleTreeGeneric::<Blake3Hash>::new(size);
            for i in 0..size {
                let mut leaf = [0u8; 32];
                leaf[0..8].copy_from_slice(&(i as u64).to_le_bytes());
                tree.set_leaf(i, leaf);
            }
            b.iter(|| {
                let mut t = tree.clone();
                black_box(t.build())
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_field_operations,
    bench_polynomial,
    bench_fft,
    bench_merkle_tree,
    bench_transcript,
    bench_ligero,
    bench_simple_proof,
    bench_larger_circuit,
    bench_hash_comparison,
    bench_merkle_hash_comparison,
);

criterion_main!(benches);
