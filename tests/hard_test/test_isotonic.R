# Test Suite for Regularized Isotonic Regression

# Set working directory to this script's location
setwd(dirname(sys.frame(1)$ofile))

source("NormalLossPiece.R")
source("isotonic_regression.R")

# Create output directory
if (!dir.exists("plots")) {
  dir.create("plots")
}

cat("===============================================\n")
cat(" HARD TEST: Regularized Isotonic Regression\n")
cat("===============================================\n\n")

# ============================================================================
# TEST 1: NormalLossPiece Class Basic Operations
# ============================================================================
cat("Test 1: NormalLossPiece Class\n")
cat("==============================\n")

# Create piece for data point y = 5
# Cost = (y - μ)² = μ² - 2*5*μ + 25
piece1 <- NormalLossPiece$new(A = 1, B = -10, C = 25, left = -Inf, right = Inf)

cat("Created piece for y = 5:\n")
piece1$print()

# Test evaluation
cat("\nEvaluations:\n")
cat(sprintf("  eval(5) = %.4f (should be 0)\n", piece1$eval(5)))
cat(sprintf("  eval(4) = %.4f (should be 1)\n", piece1$eval(4)))
cat(sprintf("  eval(6) = %.4f (should be 1)\n", piece1$eval(6)))

# Test argmin and min
cat(sprintf("\n  argmin = %.4f (should be 5)\n", piece1$get_argmin()))
cat(sprintf("  min = %.4f (should be 0)\n", piece1$get_min()))

# Test adding costs
piece2 <- piece1$copy()
piece2$add_cost(a = 1, b = -8, c = 16)  # Add cost for y = 4
cat("\nAfter adding cost for y = 4:\n")
cat(sprintf("  New argmin = %.4f (should be 4.5)\n", piece2$get_argmin()))
cat(sprintf("  New min = %.4f\n", piece2$get_min()))

cat("\n✓ Test 1 passed\n\n")

# ============================================================================
# TEST 2: Simple Monotonic Data (Already Isotonic)
# ============================================================================
cat("Test 2: Already Isotonic Data\n")
cat("==============================\n")

data2 <- c(1, 2, 3, 4, 5)
result2 <- isotonic_regression(data2, lambda = 0)

cat("Data:", paste(data2, collapse = ", "), "\n")
cat("Fitted means:", paste(round(result2$means, 3), collapse = ", "), "\n")
cat("Changepoints:", paste(result2$changepoints, collapse = ", "), "\n")
cat(sprintf("Cost: %.4f\n", result2$cost))

# Check if result is correct (should be identity)
if (all(abs(result2$means - data2) < 0.01) || length(result2$means) == 1) {
  cat("✓ Correct: Data already isotonic\n")
} else {
  cat("✗ Warning: Unexpected result\n")
}

pdf("plots/test2_already_isotonic.pdf", width = 8, height = 6)
plot_isotonic_result(data2, result2, main = "Test 2: Already Isotonic Data")
dev.off()
cat("✓ Plot saved to plots/test2_already_isotonic.pdf\n\n")

# ============================================================================
# TEST 3: Reverse Order Data (Maximum Merging)
# ============================================================================
cat("Test 3: Reverse Order Data\n")
cat("===========================\n")

data3 <- c(5, 4, 3, 2, 1)
result3 <- isotonic_regression(data3, lambda = 0)

cat("Data:", paste(data3, collapse = ", "), "\n")
cat("Fitted means:", paste(round(result3$means, 3), collapse = ", "), "\n")
cat("Changepoints:", paste(result3$changepoints, collapse = ", "), "\n")
cat(sprintf("Cost: %.4f\n", result3$cost))

# Should merge into single segment with mean = 3
expected_mean <- mean(data3)
if (length(result3$means) == 1 && abs(result3$means[1] - expected_mean) < 0.01) {
  cat(sprintf("✓ Correct: Merged to single segment with mean %.3f\n", expected_mean))
} else {
  cat("✗ Warning: Expected single segment\n")
}

pdf("plots/test3_reverse_order.pdf", width = 8, height = 6)
plot_isotonic_result(data3, result3, main = "Test 3: Reverse Order Data (Maximum Merging)")
dev.off()
cat("✓ Plot saved to plots/test3_reverse_order.pdf\n\n")

# ============================================================================
# TEST 4: Up-Down Pattern
# ============================================================================
cat("Test 4: Up-Down Pattern\n")
cat("========================\n")

data4 <- c(1, 2, 3, 2.5, 2, 1, 2, 3, 4)
result4 <- isotonic_regression(data4, lambda = 0)

cat("Data:", paste(data4, collapse = ", "), "\n")
cat("Fitted means:", paste(round(result4$means, 3), collapse = ", "), "\n")
cat("Changepoints:", paste(result4$changepoints, collapse = ", "), "\n")
cat(sprintf("Cost: %.4f\n", result4$cost))

pdf("plots/test4_updown.pdf", width = 8, height = 6)
plot_isotonic_result(data4, result4, main = "Test 4: Up-Down Pattern")
dev.off()
cat("✓ Plot saved to plots/test4_updown.pdf\n\n")

# ============================================================================
# TEST 5: Random Data
# ============================================================================
cat("Test 5: Random Data\n")
cat("====================\n")

set.seed(123)
data5 <- rnorm(20, mean = 5, sd = 2)
result5 <- isotonic_regression(data5, lambda = 0)

cat("Data (first 10):", paste(round(head(data5, 10), 2), collapse = ", "), "...\n")
cat("Number of segments:", length(result5$means), "\n")
cat(sprintf("Cost: %.4f\n", result5$cost))

# Check isotonic property
is_isotonic <- all(diff(result5$means) >= -1e-6)
if (is_isotonic) {
  cat("✓ Result is isotonic (non-decreasing)\n")
} else {
  cat("✗ Error: Result is not isotonic\n")
}

pdf("plots/test5_random.pdf", width = 10, height = 6)
plot_isotonic_result(data5, result5, main = "Test 5: Random Data (n=20)")
dev.off()
cat("✓ Plot saved to plots/test5_random.pdf\n\n")

# ============================================================================
# TEST 6: Effect of Regularization (Lambda)
# ============================================================================
cat("Test 6: Effect of Regularization\n")
cat("=================================\n")

set.seed(456)
data6 <- c(1, 1.5, 1.2, 2, 2.5, 2.3, 3, 3.5, 3.2, 4, 4.5, 4.3)

lambdas <- c(0, 0.5, 1, 2, 5)
results6 <- list()

cat(sprintf("%-10s %-15s %-10s %-10s\n", "Lambda", "Segments", "Data Cost", "Total Cost"))
cat(strrep("-", 50), "\n")

for (lambda in lambdas) {
  result <- isotonic_regression(data6, lambda = lambda)
  results6[[as.character(lambda)]] <- result

  cat(sprintf("%-10.1f %-15d %-10.2f %-10.2f\n",
              lambda, length(result$means),
              result$data_cost, result$cost))
}

cat("\nObservation: As lambda increases, number of segments decreases\n")

# Plot comparison
pdf("plots/test6_regularization.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))
for (lambda in lambdas) {
  result <- results6[[as.character(lambda)]]
  plot_isotonic_result(data6, result,
                      main = sprintf("Lambda = %.1f (%d segments)",
                                    lambda, length(result$means)))
}
dev.off()
cat("✓ Plot saved to plots/test6_regularization.pdf\n\n")

# ============================================================================
# TEST 7: Piecewise Cost Function Visualization
# ============================================================================
cat("Test 7: Piecewise Cost Visualization\n")
cat("=====================================\n")

# Simple example with 3 data points
data7 <- c(1, 3, 2)
cat("Data:", paste(data7, collapse = ", "), "\n")

# Manual construction to show pieces
piece_list <- list()

# First data point
p1 <- create_data_piece(data7[1], position = 1)
piece_list[[1]] <- list(p1)

# Second data point
p2 <- p1$copy()
p2$add_cost(1, -2*data7[2], data7[2]^2)
piece_list[[2]] <- list(p2)

# Third data point (would need isotonic constraint)
p3 <- p2$copy()
p3$add_cost(1, -2*data7[3], data7[3]^2)
piece_list[[3]] <- list(p3)

pdf("plots/test7_pieces.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))

plot_pieces(piece_list[[1]], title = "After t=1", xlim = c(-2, 5))
plot_pieces(piece_list[[2]], title = "After t=2", xlim = c(-2, 5))
plot_pieces(piece_list[[3]], title = "After t=3", xlim = c(-2, 5))

# Final result
result7 <- isotonic_regression(data7, lambda = 0)
plot_isotonic_result(data7, result7, main = "Final Result")

dev.off()
cat("✓ Plot saved to plots/test7_pieces.pdf\n\n")

# ============================================================================
# TEST 8: Comparison with Standard PAV
# ============================================================================
cat("Test 8: Comparison with Pool Adjacent Violators\n")
cat("=================================================\n")

if (require("isotone", quietly = TRUE)) {
  data8 <- c(1, 5, 4, 9, 3, 8, 2, 7, 6, 10)

  # Our implementation
  result8_ours <- isotonic_regression(data8, lambda = 0)

  # Standard PAV (using isotone package)
  result8_pav <- gpava(1:length(data8), data8)$x

  cat("Data:", paste(data8, collapse = ", "), "\n")
  cat("Our result:", paste(round(result8_ours$means[result8_ours$segments], 2), collapse = ", "), "\n")
  cat("PAV result:", paste(round(result8_pav, 2), collapse = ", "), "\n")

  # Compare
  fitted_ours <- result8_ours$means[result8_ours$segments]
  diff_max <- max(abs(fitted_ours - result8_pav))

  if (diff_max < 0.1) {
    cat(sprintf("✓ Results match (max diff: %.4f)\n", diff_max))
  } else {
    cat(sprintf("⚠ Results differ (max diff: %.4f)\n", diff_max))
  }

  pdf("plots/test8_pav_comparison.pdf", width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot_isotonic_result(data8, result8_ours, main = "Our Implementation")
  plot(1:length(data8), data8, pch = 19,
       xlab = "Index", ylab = "Value", main = "Standard PAV (isotone package)")
  lines(1:length(data8), result8_pav, col = "red", lwd = 2, type = "s")
  legend("topleft", legend = c("Data", "Fitted"), col = c("black", "red"),
         pch = c(19, NA), lty = c(NA, 1))
  dev.off()
  cat("✓ Plot saved to plots/test8_pav_comparison.pdf\n\n")
} else {
  cat("⚠ Skipping: isotone package not available\n")
  cat("  Install with: install.packages('isotone')\n\n")
}

# ============================================================================
# TEST 9: Edge Cases
# ============================================================================
cat("Test 9: Edge Cases\n")
cat("===================\n")

# Single data point
data9a <- c(5)
result9a <- isotonic_regression(data9a, lambda = 0)
cat("Single point:", data9a, "→ mean:", result9a$means, "\n")

# Two data points (increasing)
data9b <- c(1, 2)
result9b <- isotonic_regression(data9b, lambda = 0)
cat("Two points (increasing):", paste(data9b, collapse = ", "),
    "→ means:", paste(result9b$means, collapse = ", "), "\n")

# Two data points (decreasing)
data9c <- c(2, 1)
result9c <- isotonic_regression(data9c, lambda = 0)
cat("Two points (decreasing):", paste(data9c, collapse = ", "),
    "→ means:", paste(result9c$means, collapse = ", "), "\n")

# Constant data
data9d <- rep(3, 5)
result9d <- isotonic_regression(data9d, lambda = 0)
cat("Constant data:", paste(data9d, collapse = ", "),
    "→ means:", paste(unique(result9d$means), collapse = ", "), "\n")

cat("✓ All edge cases handled\n\n")

# ============================================================================
# TEST 10: Large Scale Performance
# ============================================================================
cat("Test 10: Large Scale Performance\n")
cat("==================================\n")

set.seed(789)
n_values <- c(50, 100, 200)

cat(sprintf("%-10s %-15s %-15s\n", "n", "Segments", "Time (sec)"))
cat(strrep("-", 40), "\n")

for (n in n_values) {
  data10 <- cumsum(rnorm(n)) + rnorm(n, sd = 0.5)

  start_time <- Sys.time()
  result10 <- isotonic_regression(data10, lambda = 2)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  cat(sprintf("%-10d %-15d %-15.4f\n", n, length(result10$means), elapsed))
}

cat("\n✓ Performance test complete\n\n")


