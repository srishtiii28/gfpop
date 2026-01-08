# Time-Dependent Constraints in gfpop

C++ implementations for GSoC 2025 assessment tests: https://github.com/rstats-gsoc/gsoc2025/wiki/time-dependent-constraints-in-gfpop

## Tests Completed

**Medium Test:** Unconstrained optimal partitioning with `set_to_unconstrained_min_of` method
**Hard Test:** Regularized isotonic regression with `NormalLossPiece` C++ class

---

## Medium Test: Unconstrained Optimal Partitioning

### Approach
Added `set_to_unconstrained_min_of` method to `PiecewisePoissonLossLog` class in PeakSegOptimal package:
- Copies all cost function pieces without isotonic constraints
- Enables standard changepoint detection without ordering constraints on segment means
- Simpler than `set_to_min_less_of` and `set_to_min_more_of` constrained versions

### Files

**In `PeakSegOptimal/src/`:**
- `funPieceListLog.h` - Added method declaration
- `funPieceListLog.cpp` - Method implementation
- `Makevars` - Build configuration

**Test:**
- `tests/medium_test_cpp/test_unconstrained.R` - Demonstration script

### Running the Test

```bash
cd PeakSegOptimal
R CMD INSTALL .
cd ../tests/medium_test_cpp
Rscript test_unconstrained.R
```

---

## Hard Test: NormalLossPiece C++ Class

### Approach
Implemented C++ class for Normal loss with quadratic cost functions:
- Created `NormalLossPiece` class representing A·μ² + B·μ + C
- Implemented required methods: `argmin()`, `getCost()`, `get_smaller_root()`, `get_larger_root()`
- Built isotonic regression solver using dynamic programming
- R interface via Rcpp

### Files

**In `gfpop/src/`:**
- `NormalLossPiece.h` and `NormalLossPiece.cpp` - Core class
- `IsotonicRegression.h` and `IsotonicRegression.cpp` - Solver
- `IsotonicInterface.cpp` - Rcpp interface
- `Makevars` - Build configuration

**Test:**
- `tests/hard_test_cpp/test_NormalLossPiece.R` - Demonstration script

### Running the Test

```bash
cd gfpop
R CMD INSTALL .
cd ../tests/hard_test_cpp
Rscript test_NormalLossPiece.R
```

---

## Contact

**Mentors:**
- Vincent Runge (vincent.runge@univ-evry.fr)
- Toby Dylan Hocking (tdhock5@gmail.com)
