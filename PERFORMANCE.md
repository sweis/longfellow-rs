# Performance Analysis

## Benchmark Results (longfellow-rs)

Benchmarked on: Linux 4.4.0, release build with optimizations enabled.

### Field Operations (Fp128 = 2^128 - 2^108 + 1)

| Operation | Time      | Notes                    |
|-----------|-----------|--------------------------|
| Add       | 1.58 ns   | Fast modular addition    |
| Sub       | 1.59 ns   | Fast modular subtraction |
| Mul       | 4.87 ns   | 128-bit multiplication   |
| Square    | 4.83 ns   | Optimized (3 muls vs 4)  |
| Negate    | 1.27 ns   | Simple negation          |
| Invert    | 3.54 µs   | Fermat's little theorem  |

### Batch Inversion (Montgomery's Trick) - NEW

| Elements | Batch Time | Sequential Time | Speedup |
|----------|------------|-----------------|---------|
| 10       | 4.22 µs    | 35.5 µs         | **8.4x**  |
| 50       | 4.45 µs    | 168.5 µs        | **37.8x** |
| 100      | 5.86 µs    | 341.5 µs        | **58.3x** |
| 500      | 15.6 µs    | 1.73 ms         | **111x**  |

### Polynomial Operations

| Operation              | Size | Before    | After     | Speedup   |
|------------------------|------|-----------|-----------|-----------|
| extend (coeff→eval)    | 16   | 2.70 µs   | 2.70 µs   | -         |
| extend (coeff→eval)    | 32   | 14.9 µs   | 14.9 µs   | -         |
| extend (coeff→eval)    | 64   | 64.1 µs   | 64.1 µs   | -         |
| extend (coeff→eval)    | 128  | 273 µs    | 273 µs    | -         |
| extend_evaluations     | 16   | 806 µs    | 48.7 µs   | **16.5x** |
| extend_evaluations     | 32   | 3.87 ms   | 235 µs    | **16.5x** |
| extend_evaluations     | 64   | 18.2 ms   | 1.06 ms   | **17.2x** |
| extend_evaluations     | 128  | 92.7 ms   | 5.83 ms   | **15.9x** |
| interpolate            | 8    | 168 µs    | 11.1 µs   | **15.1x** |
| interpolate            | 16   | 804 µs    | 61.1 µs   | **13.2x** |
| interpolate            | 32   | 3.83 ms   | 557 µs    | **6.9x**  |
| evaluate (32 terms)    | 32   | 151 ns    | 153 ns    | -         |

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
| 10           | 145 µs      | 69 Kelem/s       |
| 50           | 443 µs      | 113 Kelem/s      |
| 100          | 530 µs      | 189 Kelem/s      |

### ZK Proof Generation & Verification

| Circuit         | Prove Time | Verify Time |
|-----------------|------------|-------------|
| Simple (1 mul)  | 159 µs     | 19.4 ns     |
| 2 muls          | 163 µs     | 18.6 ns     |
| 5 muls          | 309 µs     | 19.1 ns     |
| 10 muls         | 842 µs     | 18.2 ns     |

## Optimizations Implemented

### 1. Batch Field Inversion (Montgomery's Trick)

**Implementation**: `field::batch_invert()`

Computes n inversions using only 1 inversion + 3(n-1) multiplications instead of n inversions.

**Impact**: Up to **111x speedup** for 500 elements. This optimization dramatically improves any operation requiring multiple inversions, particularly polynomial interpolation.

### 2. Barycentric Interpolation with Precomputed Weights

**Implementation**: `polynomial::BarycentricWeights`, `polynomial::extend_evaluations_fast()`

Uses the barycentric form of Lagrange interpolation with:
- Precomputed barycentric weights (computed once, reused many times)
- Batch inversion for the evaluation denominators
- Special case handling for evaluation at interpolation points

**Impact**: **15-17x speedup** for extend_evaluations, reducing O(n²) inversions to O(n) multiplications + O(1) inversions.

### 3. Optimized Polynomial Interpolation

**Implementation**: `Polynomial::interpolate()` now uses batch inversion

Pre-computes all (points[i] - points[j]) denominators and inverts them in a single batch operation.

**Impact**: **7-15x speedup** depending on size, with larger benefits for smaller interpolations where the overhead of individual inversions was proportionally higher.

### 4. Optimized Field Squaring

**Implementation**: `Fp128::square()`

Uses 3 multiplications instead of 4 by exploiting the symmetry in a² = a_lo² + 2*a_lo*a_hi*2^64 + a_hi²*2^128.

**Impact**: ~1% improvement (minimal due to reduction step dominating).

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

| Component          | This Library | Expected Reference | Status     |
|--------------------|--------------|-------------------|------------|
| Field Mul          | 4.87 ns      | ~3-4 ns (AVX)     | ~1.2x      |
| Field Invert       | 3.54 µs      | ~2-3 µs           | ~1.3x      |
| Batch Invert (100) | 5.86 µs      | ~5-6 µs           | Comparable |
| Polynomial Extend  | O(n²)        | O(n log n) FFT    | Needs FFT  |
| Merkle Build       | 8-9 Melem/s  | ~10-15 Melem/s    | ~0.7x      |
| Proof Gen (simple) | 159 µs       | ~100-150 µs       | ~1.1x      |

## Remaining Optimization Opportunities

### High Priority

1. **FFT for polynomial operations**
   - Would reduce `extend` from O(n²) to O(n log n)
   - Expected additional speedup: 10-100x for large polynomials
   - Requires finding primitive roots of unity for our field

2. **Parallelization**
   - Ligero commitment is embarrassingly parallel
   - Row computations can be done independently
   - Would benefit from rayon or similar

### Medium Priority

3. **SIMD for field operations**
   - AVX2/AVX-512 can process multiple field elements at once
   - Would improve all inner loops

4. **Memory allocation reduction**
   - Pre-allocate vectors where possible
   - Reuse buffers across operations

### Low Priority

5. **Montgomery form for field elements**
   - Would speed up sequences of multiplications
   - Requires conversion cost at boundaries

6. **Precomputation tables**
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

| Date       | Version | Key Improvement              | Impact               |
|------------|---------|------------------------------|----------------------|
| 2026-01-22 | 0.1.0   | Initial baseline             | -                    |
| 2026-01-22 | 0.1.1   | Batch inversion              | 111x for 500 elems   |
| 2026-01-22 | 0.1.1   | Barycentric interpolation    | 16x extend_evals     |
| 2026-01-22 | 0.1.1   | Optimized Poly::interpolate  | 15x for 8 points     |

---

*Last updated: 2026-01-22*
