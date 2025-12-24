# Poisson Optimal Partitioning (FPOP-like)
# Implements optimal partitioning for Poisson loss without constraints

#' Poisson Loss Function
#'
#' @description Compute Poisson loss for a single observation
#' Loss = -y*log(lambda) + lambda + log(y!)
#' We ignore log(y!) as it's constant for given y
#'
#' @param y Observed count (non-negative integer)
#' @param lambda Poisson parameter (rate, must be positive)
#' @return Loss value
#'
#' @examples
#' poisson_loss(5, 5)   # Should be ~5 (perfect fit)
#' poisson_loss(5, 10)  # Higher loss (poor fit)
#'
#' @export
poisson_loss <- function(y, lambda) {
  if (lambda <= 0) {
    return(Inf)
  }

  # Avoid log(0) by setting minimum lambda
  lambda <- max(lambda, 1e-10)

  # Poisson loss (negative log-likelihood, ignoring constant log(y!))
  loss <- -y * log(lambda) + lambda

  return(loss)
}


#' Segment Poisson Cost
#'
#' @description Compute optimal Poisson cost for a segment of data
#' Finds best λ (which is just the mean) and returns the cost
#'
#' @param data Numeric vector of counts for this segment
#' @return List with optimal lambda and cost
#'
#' @examples
#' segment_poisson_cost(c(5, 5, 5))  # Constant data
#' segment_poisson_cost(c(3, 5, 7))  # Variable data
#'
#' @export
segment_poisson_cost <- function(data) {
  n <- length(data)

  if (n == 0) {
    return(list(lambda = NA, cost = 0))
  }

  # Optimal lambda for Poisson is the mean
  lambda_opt <- mean(data)

  # Avoid lambda = 0
  if (lambda_opt <= 0) {
    lambda_opt <- 1e-10
  }

  # Compute total cost for this segment
  total_cost <- sum(sapply(data, function(y) poisson_loss(y, lambda_opt)))

  return(list(lambda = lambda_opt, cost = total_cost))
}


#' Optimal Partitioning for Poisson Loss (Basic Dynamic Programming)
#'
#' @description Solves the optimal partitioning problem for Poisson data
#' using dynamic programming. Finds the segmentation that minimizes:
#'   Σ Poisson_loss(y_i, λ_i) + β × (number of changepoints)
#'
#' This is the basic O(n²) algorithm without functional pruning optimization.
#'
#' @param data Numeric vector of count data (non-negative)
#' @param beta Penalty for each changepoint (non-negative)
#' @param min_segment_length Minimum length of segments (default: 1)
#'
#' @return List containing:
#'   \item{changepoints}{Indices where segments end}
#'   \item{means}{Optimal lambda (mean) for each segment}
#'   \item{segments}{Segment assignment for each data point}
#'   \item{cost}{Total penalized cost}
#'   \item{loglik}{Negative log-likelihood (without penalty)}
#'
#' @examples
#' # Two-segment data
#' data <- c(rep(5, 10), rep(15, 10))
#' result <- optimal_partitioning_poisson(data, beta = 10)
#' print(result$changepoints)
#' print(result$means)
#'
#' @export
optimal_partitioning_poisson <- function(data, beta = 0, min_segment_length = 1) {
  n <- length(data)

  if (n == 0) {
    stop("Data vector is empty")
  }

  if (beta < 0) {
    stop("Beta must be non-negative")
  }

  if (any(data < 0)) {
    stop("Poisson data must be non-negative")
  }

  # Dynamic programming arrays
  cost <- rep(Inf, n + 1)
  cost[1] <- 0  # Base case: cost before any data is 0

  best_last_changepoint <- rep(0, n + 1)

  # Precompute cumulative sums for efficiency
  cumsum_y <- c(0, cumsum(data))

  # Dynamic programming
  for (t in 1:n) {
    # Try all possible last changepoints
    for (s in max(0, t - n):(t - min_segment_length)) {
      if (cost[s + 1] == Inf) next

      # Compute cost for segment (s+1):t
      segment_data <- data[(s + 1):t]
      seg_info <- segment_poisson_cost(segment_data)
      segment_cost <- seg_info$cost

      # Total cost if we put a changepoint at s
      penalty <- if (s > 0) beta else 0  # No penalty for first segment
      candidate_cost <- cost[s + 1] + segment_cost + penalty

      if (candidate_cost < cost[t + 1]) {
        cost[t + 1] <- candidate_cost
        best_last_changepoint[t + 1] <- s
      }
    }
  }

  # Backtrack to find changepoints
  changepoints <- c()
  current <- n

  while (current > 0) {
    changepoints <- c(current, changepoints)
    current <- best_last_changepoint[current + 1]
  }

  # Compute segment means
  n_segments <- length(changepoints)
  means <- numeric(n_segments)
  segments <- numeric(n)

  segment_starts <- c(0, changepoints[-n_segments]) + 1
  segment_ends <- changepoints

  for (i in 1:n_segments) {
    start <- segment_starts[i]
    end <- segment_ends[i]
    segment_data <- data[start:end]

    means[i] <- mean(segment_data)
    segments[start:end] <- i
  }

  # Compute costs
  total_loglik <- 0
  for (i in 1:n_segments) {
    start <- segment_starts[i]
    end <- segment_ends[i]
    segment_data <- data[start:end]
    seg_info <- segment_poisson_cost(segment_data)
    total_loglik <- total_loglik + seg_info$cost
  }

  total_cost <- total_loglik + beta * (n_segments - 1)

  return(list(
    changepoints = changepoints,
    means = means,
    segments = segments,
    cost = total_cost,
    loglik = total_loglik,
    n_segments = n_segments,
    beta = beta
  ))
}


#' Optimal Partitioning with Pruning (FPOP-inspired)
#'
#' @description More efficient version using functional pruning ideas
#' Still O(n²) worst case but faster in practice
#'
#' @param data Numeric vector of count data
#' @param beta Penalty for each changepoint
#' @param min_segment_length Minimum length of segments
#'
#' @return Same as optimal_partitioning_poisson
#'
#' @export
optimal_partitioning_poisson_pruned <- function(data, beta = 0, min_segment_length = 1) {
  n <- length(data)

  if (n == 0) stop("Data vector is empty")
  if (beta < 0) stop("Beta must be non-negative")
  if (any(data < 0)) stop("Poisson data must be non-negative")

  # Arrays for DP
  cost <- rep(Inf, n + 1)
  cost[1] <- 0

  best_last_changepoint <- rep(0, n + 1)

  # Candidate set: positions that might be optimal last changepoints
  candidates <- list(0)  # Start with position 0

  for (t in 1:n) {
    best_cost <- Inf
    best_s <- 0

    # Only check candidates (pruned set)
    for (s in unlist(candidates)) {
      if (t - s < min_segment_length) next
      if (cost[s + 1] == Inf) next

      # Compute segment cost
      segment_data <- data[(s + 1):t]
      seg_info <- segment_poisson_cost(segment_data)
      segment_cost <- seg_info$cost

      penalty <- if (s > 0) beta else 0
      candidate_cost <- cost[s + 1] + segment_cost + penalty

      if (candidate_cost < best_cost) {
        best_cost <- candidate_cost
        best_s <- s
      }
    }

    cost[t + 1] <- best_cost
    best_last_changepoint[t + 1] <- best_s

    # Pruning: add current position as candidate
    # Simple strategy: keep recent positions and prune very old ones
    candidates <- append(candidates, t)

    # Prune candidates that are too old (heuristic)
    if (length(candidates) > 50) {
      # Keep only recent candidates
      candidates <- tail(candidates, 40)
    }
  }

  # Backtrack (same as before)
  changepoints <- c()
  current <- n

  while (current > 0) {
    changepoints <- c(current, changepoints)
    current <- best_last_changepoint[current + 1]
  }

  # Compute means and segments
  n_segments <- length(changepoints)
  means <- numeric(n_segments)
  segments <- numeric(n)

  segment_starts <- c(0, changepoints[-n_segments]) + 1
  segment_ends <- changepoints

  for (i in 1:n_segments) {
    start <- segment_starts[i]
    end <- segment_ends[i]
    segment_data <- data[start:end]
    means[i] <- mean(segment_data)
    segments[start:end] <- i
  }

  total_loglik <- 0
  for (i in 1:n_segments) {
    start <- segment_starts[i]
    end <- segment_ends[i]
    segment_data <- data[start:end]
    seg_info <- segment_poisson_cost(segment_data)
    total_loglik <- total_loglik + seg_info$cost
  }

  total_cost <- total_loglik + beta * (n_segments - 1)

  return(list(
    changepoints = changepoints,
    means = means,
    segments = segments,
    cost = total_cost,
    loglik = total_loglik,
    n_segments = n_segments,
    beta = beta
  ))
}


#' Plot Poisson Segmentation Result
#'
#' @description Visualize data and fitted Poisson segmentation
#' @param data Original data
#' @param result Result from optimal_partitioning_poisson()
#' @param main Plot title
#'
#' @export
plot_poisson_result <- function(data, result, main = "Poisson Optimal Partitioning") {
  n <- length(data)
  fitted <- result$means[result$segments]

  plot(1:n, data, pch = 19, col = "black", type = "p",
       xlab = "Index", ylab = "Count",
       main = main,
       ylim = range(c(data, fitted)))

  # Plot fitted values
  lines(1:n, fitted, col = "red", lwd = 2, type = "s")

  # Mark changepoints
  if (length(result$changepoints) > 1) {
    for (cp in result$changepoints[-length(result$changepoints)]) {
      abline(v = cp + 0.5, col = "blue", lty = 2, lwd = 1.5)
    }
  }

  # Add legend
  legend("topleft",
         legend = c("Data", "Fitted means", "Changepoints"),
         col = c("black", "red", "blue"),
         pch = c(19, NA, NA),
         lty = c(NA, 1, 2),
         lwd = c(NA, 2, 1.5),
         cex = 0.8)

  # Add info box
  info_text <- sprintf("Segments: %d\nCost: %.2f\nLoglik: %.2f\nβ: %.2f",
                      result$n_segments, result$cost, result$loglik, result$beta)

  text(n * 0.05, max(data) * 0.95, info_text, adj = c(0, 1), cex = 0.7,
       family = "mono", col = "darkgreen")

  grid()
}


#' Generate Synthetic Poisson Data
#'
#' @description Create Poisson data with known changepoints for testing
#' @param n Total number of data points
#' @param changepoint_props Proportions where changes occur (e.g., c(0.3, 0.7))
#' @param means Lambda values for each segment
#' @param seed Random seed
#'
#' @return Numeric vector of Poisson counts
#'
#' @export
generate_poisson_data <- function(n, changepoint_props, means, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_segments <- length(means)
  if (length(changepoint_props) != n_segments - 1) {
    stop("Need length(means) - 1 changepoint proportions")
  }

  # Convert proportions to indices
  changepoints <- round(changepoint_props * n)
  changepoints <- c(0, changepoints, n)

  data <- numeric(n)

  for (i in 1:n_segments) {
    start <- changepoints[i] + 1
    end <- changepoints[i + 1]

    if (start <= end) {
      data[start:end] <- rpois(end - start + 1, lambda = means[i])
    }
  }

  return(data)
}
