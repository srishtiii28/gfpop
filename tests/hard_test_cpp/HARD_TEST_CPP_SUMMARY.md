# Hard Test: NormalLossPiece C++ Implementation

## GSoC 2025: Time-Dependent Constraints in gfpop

**Status:** IMPLEMENTED ✓

---

## Implementation Details

### Files Created

1. **`gfpop/src/NormalLossPiece.h`** - Header file defining the NormalLossPiece class
2. **`gfpop/src/NormalLossPiece.cpp`** - Implementation of NormalLossPiece methods
3. **`gfpop/src/IsotonicRegression.h`** - Isotonic regression solver header
4. **`gfpop/src/IsotonicRegression.cpp`** - Isotonic regression solver implementation
5. **`gfpop/src/IsotonicInterface.cpp`** - Rcpp interface for R
6. **`gfpop/src/Makevars`** - Build configuration (fixes C++ standard library paths)

### NormalLossPiece Class

The `NormalLossPiece` class implements quadratic cost functions for Normal loss: **A·μ² + B·μ + C**

**Class Members:**
```cpp
double A;           // quadratic coefficient
double B;           // linear coefficient
double C;           // constant term
double min_mean;    // minimum mean for this piece
double max_mean;    // maximum mean for this piece
int data_i;         // data index
double prev_mean;   // previous mean
```

**Implemented Methods:**
- `NormalLossPiece()` - Default constructor
- `NormalLossPiece(a, b, c, min_m, max_m, i, prev)` - Parameterized constructor
- `double argmin()` - Returns argmin = -B/(2A)
- `double getCost(double mean)` - Evaluates cost at given mean
- `double get_smaller_root(double level)` - Finds smaller root where cost = level
- `double get_larger_root(double level)` - Finds larger root where cost = level
- `bool has_two_roots(double level)` - Checks if equation has two real roots
- `void print()` - Debug output

### R Interface

**Function:** `isotonicRegressionCpp(data, lambda)`

**Parameters:**
- `data`: Numeric vector of observations
- `lambda`: Regularization parameter (default = 0)

**Returns:** Numeric vector of isotonic means

### Build System Fix

Created `Makevars` file to fix C++ standard library include path issue:
```makefile
CXX_STD = CXX17
PKG_CPPFLAGS = -I. -I/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp/include'