# Performance Analysis

## Benchmark Results (longfellow-rs)

Benchmarked on: Linux 4.4.0, release build with optimizations enabled.

### Field Operations (Fp128 = 2^128 - 2^108 + 1)

| Operation | Time      | Notes                    |
|-----------|-----------|--------------------------|
| Add       | 1.58 ns   | Fast modular addition    |
| Sub       | 1.59 ns   | Fast modular subtraction |
| Mul       | 4.87 ns   | 128-bit multiplication   |
| Square    | 4.90 ns   | Same as mul (not optimized) |
| Negate    | 1.27 ns   | Simple negation          |
| Invert    | 3.54 µs   | Fermat's little theorem  |

### Polynomial Operations

| Operation              | Size | Time      | Throughput      |
|------------------------|------|-----------|-----------------|
| extend (coeff→eval)    | 16   | 2.70 µs   | 5.93 Melem/s    |
| extend (coeff→eval)    | 32   | 14.9 µs   | 2.14 Melem/s    |
| extend (coeff→eval)    | 64   | 64.1 µs   | 999 Kelem/s     |
| extend (coeff→eval)    | 128  | 273 µs    | 470 Kelem/s     |
| extend_evaluations     | 16   | 806 µs    | 19.9 Kelem/s    |
| extend_evaluations     | 32   | 3.87 ms   | 8.27 Kelem/s    |
| extend_evaluations     | 64   | 18.2 ms   | 3.52 Kelem/s    |
| extend_evaluations     | 128  | 92.7 ms   | 1.38 Kelem/s    |
| interpolate            | 8    | 168 µs    | 752 Kelem/s     |
| interpolate            | 16   | 804 µs    | 159 Kelem/s     |
| interpolate            | 32   | 3.83 ms   | 33.4 Kelem/s    |
| evaluate (32 terms)    | 32   | 151 ns    | 847 Melem/s     |

### Merkle Tree (SHA-256)

| Operation | Leaves | Time      | Throughput      |
|-----------|--------|-----------|-----------------|
| build     | 16     | 1.78 µs   | 9.01 Melem/s    |
| build     | 64     | 7.43 µs   | 8.61 Melem/s    |
| build     | 256    | 30.0 µs   | 8.52 Melem/s    |
| build     | 1024   | 123 µs    | 8.34 Melem/s    |
| proof     | 16     | 76.3 ns   | 210 Melem/s     |
| proof     | 64     | 176 ns    | 364 Melem/s     |
| proof     | 256    | 516 ns    | 496 Melem/s     |
| proof     | 1024   | 1.80 µs   | 570 Melem/s     |

### Transcript (Fiat-Shamir)

| Operation              | Time      |
|------------------------|-----------|
| generate_challenge     | 82.0 ns   |
| generate_challenges(10)| 819 ns    |
| write_field_element    | 60.0 ns   |

### Ligero Commitment

| Witness Size | Commit Time | Throughput       |
|--------------|-------------|------------------|
| 10           | 133 µs      | 75.2 Kelem/s     |
| 50           | 402 µs      | 124 Kelem/s      |
| 100          | 490 µs      | 204 Kelem/s      |

### ZK Proof Generation & Verification

| Circuit         | Prove Time | Verify Time |
|-----------------|------------|-------------|
| Simple (1 mul)  | 154 µs     | 17.0 ns     |
| 2 muls          | 159 µs     | 16.5 ns     |
| 5 muls          | 284 µs     | 17.0 ns     |
| 10 muls         | 798 µs     | 16.6 ns     |

## Comparison with Reference Implementations

### Reference Implementations

1. **google/longfellow-zk** (C++)
   - Reference implementation
   - Optimized for x86-64 with AVX2/AVX-512 intrinsics
   - Uses OpenMP for parallelization

2. **abetterinternet/zk-cred-longfellow** (Rust)
   - Production implementation by ISRG
   - Used in Prio system for privacy-preserving metrics
   - Optimized for batch operations

### Expected Performance Comparison

| Component          | This Library | Expected Reference | Ratio |
|--------------------|--------------|-------------------|-------|
| Field Mul          | 4.87 ns      | ~3-4 ns (AVX)     | ~1.2x |
| Field Invert       | 3.54 µs      | ~2-3 µs           | ~1.3x |
| Polynomial Extend  | O(n²)        | O(n log n) FFT    | varies |
| Merkle Build       | 8-9 Melem/s  | ~10-15 Melem/s    | ~0.7x |
| Proof Gen (simple) | 154 µs       | ~100-150 µs       | ~1.0x |

### Key Observations

1. **Field arithmetic** is competitive with reference implementations
2. **Polynomial operations** use O(n²) naive algorithms vs O(n log n) FFT
3. **extend_evaluations** is a major bottleneck due to O(n²) interpolation
4. **Merkle tree** performance is reasonable
5. **Verification** is very fast (< 20 ns for simple circuits)

## Bottleneck Analysis

### Top Performance Bottlenecks

1. **extend_evaluations** - O(n²) Lagrange interpolation
   - 128 points: 92.7 ms
   - This is 340x slower than extend() for the same size
   - Solution: Use FFT-based interpolation

2. **Polynomial interpolation** - O(n²) naive algorithm
   - 32 points: 3.83 ms
   - Solution: Use FFT with precomputed roots of unity

3. **Field inversion** - 3.54 µs per inversion
   - Uses Fermat's little theorem: a^(p-2) mod p
   - 127 squarings + 126 multiplications
   - Solution: Batch inversions using Montgomery's trick

## Optimization Opportunities

### High Priority

1. **FFT for polynomial operations**
   - Would reduce extend_evaluations from O(n²) to O(n log n)
   - Expected speedup: 10-100x for large polynomials

2. **Batch field inversions**
   - Montgomery's trick: n inversions in 1 inversion + 3(n-1) muls
   - Would help polynomial interpolation significantly

3. **Parallelization**
   - Ligero commitment is embarrassingly parallel
   - Row computations can be done independently

### Medium Priority

4. **SIMD for field operations**
   - AVX2/AVX-512 can process multiple field elements at once
   - Would improve all inner loops

5. **Squaring optimization**
   - Currently square() = mul(x, x)
   - Dedicated squaring can be faster

6. **Memory allocation reduction**
   - Pre-allocate vectors where possible
   - Reuse buffers across operations

### Low Priority

7. **Montgomery form for field elements**
   - Would speed up sequences of multiplications
   - Requires conversion cost at boundaries

8. **Precomputation tables**
   - For frequently used values in polynomial evaluation

## Running Benchmarks

```bash
# Run all benchmarks
cargo bench

# Run specific benchmark group
cargo bench -- "Field"
cargo bench -- "Polynomial"
cargo bench -- "Simple Proof"

# Generate HTML report
cargo bench -- --save-baseline current
```

## Tracking Progress

| Date       | Version | Simple Proof (prove) | Notes              |
|------------|---------|---------------------|---------------------|
| 2026-01-22 | 0.1.0   | 154 µs              | Initial baseline    |

---

*Last updated: 2026-01-22*
