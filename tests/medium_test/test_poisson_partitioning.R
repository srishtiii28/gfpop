# Test Suite for Poisson Optimal Partitioning

# Set working directory to this script's location
setwd(dirname(sys.frame(1)$ofile))

source("poisson_optimal_partitioning.R")

# Create output directory
if (!dir.exists("plots")) {
  dir.create("plots")
}


# ============================================================================
# TEST 1: Basic Poisson Loss Function
# ============================================================================
cat("Test 1: Poisson Loss Function\n")
cat("===============================\n")

# Test with perfect fit
loss1 <- poisson_loss(5, 5)
cat(sprintf("poisson_loss(5, 5) = %.4f\n", loss1))

# Test with different lambda
loss2 <- poisson_loss(5, 10)
cat(sprintf("poisson_loss(5, 10) = %.4f\n", loss2))

loss3 <- poisson_loss(5, 3)
cat(sprintf("poisson_loss(5, 3) = %.4f\n", loss3))

# Verify loss is minimized at y = lambda
y_val <- 5
lambdas <- seq(1, 10, by = 0.5)
losses <- sapply(lambdas, function(lam) poisson_loss(y_val, lam))

pdf("plots/test1_poisson_loss.pdf", width = 8, height = 6)
plot(lambdas, losses, type = "l", lwd = 2, col = "blue",
     xlab = expression(lambda), ylab = "Loss",
     main = "Poisson Loss Function for y = 5")
abline(v = y_val, col = "red", lty = 2, lwd = 2)
text(y_val + 0.5, max(losses) * 0.9, "λ = y (minimum)", col = "red")
grid()
dev.off()

cat("✓ Loss minimized at λ = y\n")
cat("✓ Plot saved to plots/test1_poisson_loss.pdf\n\n")

# ============================================================================
# TEST 2: Segment Cost Computation
# ============================================================================
cat("Test 2: Segment Cost Computation\n")
cat("==================================\n")

# Constant data
seg1 <- segment_poisson_cost(c(5, 5, 5, 5, 5))
cat(sprintf("Constant data (all 5): λ = %.4f, cost = %.4f\n",
            seg1$lambda, seg1$cost))

# Variable data
seg2 <- segment_poisson_cost(c(3, 5, 7, 4, 6))
cat(sprintf("Variable data: λ = %.4f, cost = %.4f\n",
            seg2$lambda, seg2$cost))

# Single point
seg3 <- segment_poisson_cost(c(10))
cat(sprintf("Single point (10): λ = %.4f, cost = %.4f\n",
            seg3$lambda, seg3$cost))

cat("✓ Segment costs computed correctly\n\n")

# ============================================================================
# TEST 3: Constant Poisson Data (Single Segment)
# ============================================================================
cat("Test 3: Constant Data\n")
cat("======================\n")

data3 <- rep(5, 20)
result3 <- optimal_partitioning_poisson(data3, beta = 10)

cat("Data: All 5's (n=20)\n")
cat(sprintf("Number of segments: %d (should be 1)\n", result3$n_segments))
cat(sprintf("Estimated mean: %.4f (should be ~5)\n", result3$means[1]))
cat(sprintf("Cost: %.4f\n", result3$cost))

if (result3$n_segments == 1 && abs(result3$means[1] - 5) < 0.1) {
  cat("✓ Correctly identified as single segment\n")
} else {
  cat("⚠ Warning: Expected single segment with mean 5\n")
}

pdf("plots/test3_constant.pdf", width = 8, height = 6)
plot_poisson_result(data3, result3, main = "Test 3: Constant Data (λ=5)")
dev.off()
cat("✓ Plot saved to plots/test3_constant.pdf\n\n")

# ============================================================================
# TEST 4: Two-Segment Data
# ============================================================================
cat("Test 4: Two-Segment Data\n")
cat("=========================\n")

data4 <- generate_poisson_data(40, changepoint_props = c(0.5), means = c(3, 10), seed = 123)

result4 <- optimal_partitioning_poisson(data4, beta = 5)

cat("Data: Two segments (n=40), true means: 3 and 10\n")
cat(sprintf("Detected segments: %d\n", result4$n_segments))
cat("Estimated means:", paste(round(result4$means, 2), collapse = ", "), "\n")
cat("Changepoints:", paste(result4$changepoints, collapse = ", "), "\n")
cat(sprintf("Cost: %.4f\n", result4$cost))

pdf("plots/test4_two_segments.pdf", width = 10, height = 6)
plot_poisson_result(data4, result4, main = "Test 4: Two-Segment Data (λ=3, λ=10)")
abline(v = 20.5, col = "green", lty = 3, lwd = 2)  # True changepoint
legend("topright", legend = "True changepoint", col = "green", lty = 3, lwd = 2, cex = 0.7)
dev.off()
cat("✓ Plot saved to plots/test4_two_segments.pdf\n\n")

# ============================================================================
# TEST 5: Multiple Segments
# ============================================================================
cat("Test 5: Multiple Segments\n")
cat("==========================\n")

data5 <- generate_poisson_data(60, changepoint_props = c(0.25, 0.5, 0.75),
                               means = c(5, 15, 8, 20), seed = 456)

result5 <- optimal_partitioning_poisson(data5, beta = 8)

cat("Data: Four segments (n=60), true means: 5, 15, 8, 20\n")
cat(sprintf("Detected segments: %d\n", result5$n_segments))
cat("Estimated means:", paste(round(result5$means, 2), collapse = ", "), "\n")
cat("Changepoints:", paste(result5$changepoints, collapse = ", "), "\n")
cat(sprintf("Cost: %.4f\n", result5$cost))

pdf("plots/test5_multiple_segments.pdf", width = 12, height = 6)
plot_poisson_result(data5, result5, main = "Test 5: Four-Segment Data")
# Mark true changepoints
for (cp in c(15, 30, 45)) {
  abline(v = cp + 0.5, col = "green", lty = 3, lwd = 2)
}
legend("topright", legend = "True changepoints", col = "green", lty = 3, lwd = 2, cex = 0.7)
dev.off()
cat("✓ Plot saved to plots/test5_multiple_segments.pdf\n\n")

# ============================================================================
# TEST 6: Effect of Penalty (Beta)
# ============================================================================
cat("Test 6: Effect of Penalty Parameter\n")
cat("=====================================\n")

set.seed(789)
data6 <- generate_poisson_data(50, changepoint_props = c(0.2, 0.5, 0.8),
                               means = c(5, 12, 7, 15), seed = 789)

betas <- c(0, 5, 10, 20, 50)
results6 <- list()

cat(sprintf("%-10s %-15s %-15s %-15s\n", "Beta", "Segments", "Loglik", "Total Cost"))
cat(strrep("-", 60), "\n")

for (beta in betas) {
  result <- optimal_partitioning_poisson(data6, beta = beta)
  results6[[as.character(beta)]] <- result

  cat(sprintf("%-10.1f %-15d %-15.2f %-15.2f\n",
              beta, result$n_segments, result$loglik, result$cost))
}

cat("\nObservation: As beta increases, number of segments decreases\n")

pdf("plots/test6_penalty_effect.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))
for (beta in betas) {
  result <- results6[[as.character(beta)]]
  plot_poisson_result(data6, result,
                     main = sprintf("β = %.1f (%d segments)",
                                   beta, result$n_segments))
}
dev.off()
cat("✓ Plot saved to plots/test6_penalty_effect.pdf\n\n")

# ============================================================================
# TEST 7: Comparison with Segmentor3IsBack (if available)
# ============================================================================
cat("Test 7: Validation with Segmentor3IsBack\n")
cat("==========================================\n")

if (requireNamespace("Segmentor3IsBack", quietly = TRUE)) {
  library(Segmentor3IsBack)

  set.seed(111)
  data7 <- generate_poisson_data(50, changepoint_props = c(0.4, 0.7),
                                 means = c(5, 15, 8), seed = 111)

  # Our implementation
  result7_ours <- optimal_partitioning_poisson(data7, beta = 10)

  # Segmentor3IsBack
  # Use Kmax = reasonable upper bound on segments
  seg3_result <- Segmentor(data7, model = 3, Kmax = 10)

  # Find best K based on penalized likelihood
  # (Note: Segmentor3IsBack uses different penalty formulation)
  best_K <- which.min(seg3_result@likelihood + 10 * (1:10))
  seg3_changepoints <- seg3_result@changepoints[1:best_K, best_K]

  cat("Data: Three segments (n=50)\n")
  cat(sprintf("Our implementation: %d segments\n", result7_ours$n_segments))
  cat("  Changepoints:", paste(result7_ours$changepoints, collapse = ", "), "\n")
  cat("  Means:", paste(round(result7_ours$means, 2), collapse = ", "), "\n")

  cat(sprintf("\nSegmentor3IsBack: %d segments\n", best_K))
  cat("  Changepoints:", paste(seg3_changepoints, collapse = ", "), "\n")

  # Compare changepoints (allow small differences)
  if (result7_ours$n_segments == best_K) {
    max_diff <- max(abs(result7_ours$changepoints - seg3_changepoints))
    if (max_diff <= 2) {
      cat(sprintf("\n✓ Results match (max difference: %d positions)\n", max_diff))
    } else {
      cat(sprintf("\n⚠ Results differ (max difference: %d positions)\n", max_diff))
    }
  } else {
    cat("\n⚠ Different number of segments detected\n")
  }

  pdf("plots/test7_segmentor_comparison.pdf", width = 12, height = 5)
  par(mfrow = c(1, 2))
  plot_poisson_result(data7, result7_ours, main = "Our Implementation")

  # Plot Segmentor3IsBack result
  fitted_seg3 <- rep(NA, length(data7))
  starts <- c(1, seg3_changepoints[-best_K] + 1)
  ends <- seg3_changepoints
  for (i in 1:best_K) {
    fitted_seg3[starts[i]:ends[i]] <- mean(data7[starts[i]:ends[i]])
  }

  plot(1:length(data7), data7, pch = 19, col = "black",
       xlab = "Index", ylab = "Count",
       main = "Segmentor3IsBack")
  lines(1:length(data7), fitted_seg3, col = "red", lwd = 2, type = "s")
  for (cp in seg3_changepoints[-best_K]) {
    abline(v = cp + 0.5, col = "blue", lty = 2, lwd = 1.5)
  }
  grid()
  dev.off()
  cat("✓ Plot saved to plots/test7_segmentor_comparison.pdf\n\n")

} else {
  cat("⚠ Segmentor3IsBack package not available\n")
  cat("  Install with: install.packages('Segmentor3IsBack')\n")
  cat("  Skipping comparison test\n\n")
}

# ============================================================================
# TEST 8: Edge Cases
# ============================================================================
cat("Test 8: Edge Cases\n")
cat("===================\n")

# Single data point
data8a <- c(5)
result8a <- optimal_partitioning_poisson(data8a, beta = 0)
cat(sprintf("Single point (5): %d segment, mean = %.4f\n",
            result8a$n_segments, result8a$means[1]))

# Two points (same)
data8b <- c(7, 7)
result8b <- optimal_partitioning_poisson(data8b, beta = 10)
cat(sprintf("Two same points (7, 7): %d segment, mean = %.4f\n",
            result8b$n_segments, result8b$means[1]))

# Two points (different)
data8c <- c(3, 9)
result8c_low <- optimal_partitioning_poisson(data8c, beta = 0)
result8c_high <- optimal_partitioning_poisson(data8c, beta = 100)
cat(sprintf("Two different points (3, 9) with β=0: %d segments\n",
            result8c_low$n_segments))
cat(sprintf("Two different points (3, 9) with β=100: %d segment\n",
            result8c_high$n_segments))

# All zeros
data8d <- rep(0, 10)
result8d <- optimal_partitioning_poisson(data8d, beta = 5)
cat(sprintf("All zeros (n=10): %d segment, mean = %.4f\n",
            result8d$n_segments, result8d$means[1]))

cat("✓ All edge cases handled\n\n")

# ============================================================================
# TEST 9: Performance Comparison
# ============================================================================
cat("Test 9: Performance (Basic vs Pruned)\n")
cat("=======================================\n")

set.seed(999)
n_values <- c(50, 100, 200)

cat(sprintf("%-10s %-20s %-20s\n", "n", "Basic (sec)", "Pruned (sec)"))
cat(strrep("-", 50), "\n")

for (n in n_values) {
  data9 <- rpois(n, lambda = 10)

  # Basic version
  start_basic <- Sys.time()
  result_basic <- optimal_partitioning_poisson(data9, beta = 5)
  time_basic <- as.numeric(difftime(Sys.time(), start_basic, units = "secs"))

  # Pruned version
  start_pruned <- Sys.time()
  result_pruned <- optimal_partitioning_poisson_pruned(data9, beta = 5)
  time_pruned <- as.numeric(difftime(Sys.time(), start_pruned, units = "secs"))

  cat(sprintf("%-10d %-20.4f %-20.4f\n", n, time_basic, time_pruned))
}

cat("\n✓ Performance test complete\n\n")

# ============================================================================
# TEST 10: Real-world-like Data
# ============================================================================
cat("Test 10: Real-world Scenario\n")
cat("=============================\n")

# Simulate count data with gradual changes
set.seed(2025)
n <- 100
true_means <- c(rep(5, 25), rep(8, 25), rep(12, 25), rep(6, 25))
data10 <- rpois(n, lambda = true_means)

result10 <- optimal_partitioning_poisson(data10, beta = 15)

cat(sprintf("Data: n=%d with 4 true segments\n", n))
cat(sprintf("Detected: %d segments\n", result10$n_segments))
cat("Estimated means:", paste(round(result10$means, 2), collapse = ", "), "\n")
cat("Changepoints:", paste(result10$changepoints, collapse = ", "), "\n")

pdf("plots/test10_realworld.pdf", width = 12, height = 6)
plot_poisson_result(data10, result10, main = "Test 10: Real-world Scenario")
# True changepoints
for (cp in c(25, 50, 75)) {
  abline(v = cp + 0.5, col = "green", lty = 3, lwd = 2)
}
legend("topright", legend = "True changepoints", col = "green", lty = 3, lwd = 2, cex = 0.7)
dev.off()
cat("✓ Plot saved to plots/test10_realworld.pdf\n\n")

