Assessment test implementations for the GSoC 2025 project: https://github.com/rstats-gsoc/gsoc2025/wiki/time-dependent-constraints-in-gfpop

## Tests Completed

**Medium Test:** Poisson optimal partitioning
**Hard Test:** Regularized isotonic regression

---

## Medium Test: Poisson Optimal Partitioning

### Approach
Implemented a dynamic programming solver for optimal partitioning with Poisson loss:
- Basic O(n²) algorithm with penalty-based regularization
- FPOP-inspired pruning to improve efficiency
- Validated against Segmentor3IsBack package

### Running the Test
```bash
cd tests/medium_test
Rscript test_poisson_partitioning.R
```

**Prerequisites:** `Segmentor3IsBack` package for validation

---

## Hard Test: Regularized Isotonic Regression

### Approach
Implemented an isotonic regression solver with regularization:
- Created `NormalLossPiece` R6 class for quadratic cost functions (A·μ² + B·μ + C)
- Dynamic programming with isotonic constraints (non-decreasing means)
- Min-envelope pruning for optimization
- Validated against Pool Adjacent Violators (PAV) algorithm from `isotone` package

### Running the Test
```bash
cd tests/hard_test
Rscript test_isotonic.R
```

**Prerequisites:** `R6`, `isotone` packages
