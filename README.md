# GSoC 2025: Time-Dependent Constraints in gfpop

Assessment test implementation for the GSoC 2025 project: https://github.com/rstats-gsoc/gsoc2025/wiki/time-dependent-constraints-in-gfpop

## Test Completed

**Hard Test:** Regularized isotonic regression with NormalLossPiece C++ class

---

## Hard Test: NormalLossPiece C++ Implementation

### Approach
Implemented a C++ class for Normal loss with quadratic cost functions:
- Created `NormalLossPiece` class representing A·μ² + B·μ + C
- Implemented required methods: `argmin()`, `getCost()`, `get_smaller_root()`, `get_larger_root()`
- Built isotonic regression solver using dynamic programming
- R interface via Rcpp for testing

### Files

**In `gfpop/src/`:**
- `NormalLossPiece.h` and `NormalLossPiece.cpp` - Core class implementation
- `IsotonicRegression.h` and `IsotonicRegression.cpp` - Solver implementation
- `IsotonicInterface.cpp` - Rcpp interface
- `Makevars` - Build configuration

**Tests:**
- `tests/hard_test_cpp/test_NormalLossPiece.R` - Test script
- `tests/hard_test_cpp/HARD_TEST_CPP_SUMMARY.md` - Implementation details

### Running the Test

```bash
# Build and install the modified gfpop package
cd gfpop
R CMD INSTALL .

# Run tests
cd ../tests/hard_test_cpp
Rscript test_NormalLossPiece.R
```

**Prerequisites:** `Rcpp`, `R6` packages

---

## Contact

**Mentors:**
- Vincent Runge (vincent.runge@univ-evry.fr)
- Toby Dylan Hocking (tdhock5@gmail.com)
