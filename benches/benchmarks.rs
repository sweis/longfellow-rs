//! Benchmarks for the Longfellow ZK library.

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use longfellow_zk::{
    circuit::{CircuitBuilder, LayerBuilder},
    field::Fp128,
    zk::{ZkProver, verify_zk},
    transcript::Transcript,
    merkle::MerkleTree,
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

fn bench_field_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("Field Operations");

    let a = Fp128::from_u64(12345678901234567890);
    let b = Fp128::from_u64(98765432109876543210);

    group.bench_function("add", |bencher| {
        bencher.iter(|| black_box(a) + black_box(b))
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

    group.finish();
}

fn bench_merkle_tree(c: &mut Criterion) {
    let mut group = c.benchmark_group("Merkle Tree");

    for size in [16, 64, 256, 1024].iter() {
        group.bench_with_input(BenchmarkId::new("build", size), size, |b, &size| {
            let mut tree = MerkleTree::new(size);
            for i in 0..size {
                tree.set_leaf(i, [i as u8; 32]);
            }
            b.iter(|| {
                let mut t = tree.clone();
                black_box(t.build())
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

    group.bench_function("write_and_challenge", |bencher| {
        bencher.iter(|| {
            let mut transcript = Transcript::new(b"benchmark");
            transcript.write_bytes(&[0u8; 64]);
            let c: Fp128 = transcript.generate_challenge();
            black_box(c)
        })
    });

    group.finish();
}

fn bench_simple_proof(c: &mut Criterion) {
    let mut group = c.benchmark_group("Simple Proof");

    // Create a simple multiplication circuit
    let mut builder = CircuitBuilder::<Fp128>::new(0, 2);
    let mut layer = LayerBuilder::new(0, 1);
    layer.add_mul(0, 0, 1, Fp128::ONE);
    builder.add_layer(layer);
    let circuit = builder.build();

    let witness = vec![Fp128::from_u64(3), Fp128::from_u64(4)];

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

criterion_group!(
    benches,
    bench_field_operations,
    bench_merkle_tree,
    bench_transcript,
    bench_simple_proof,
);

criterion_main!(benches);
