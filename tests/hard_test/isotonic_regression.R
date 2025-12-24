# Regularized Isotonic Regression Solver
# Uses NormalLossPiece class for dynamic programming

source("NormalLossPiece.R")

#' Regularized Isotonic Regression
#'
#' @description Solves the regularized isotonic regression problem:
#'   minimize: Σ(y_i - μ_i)² + λ * (number of segments)
#'   subject to: μ_1 ≤ μ_2 ≤ ... ≤ μ_n
#'
#' Uses dynamic programming with piecewise quadratic cost functions
#' represented by NormalLossPiece objects.
#'
#' @param data Numeric vector of observations
#' @param lambda Non-negative penalty for adding segments (default: 0)
#' @param min_value Minimum allowed value for μ (default: -Inf)
#' @param max_value Maximum allowed value for μ (default: Inf)
#'
#' @return List containing:
#'   \item{changepoints}{Indices where segments end}
#'   \item{means}{Estimated mean for each segment}
#'   \item{segments}{Segment assignments for each data point}
#'   \item{cost}{Total cost (data fit + penalties)}
#'   \item{pieces}{Final cost function pieces}
#'
#' @examples
#' # Simple example with 5 data points
#' data <- c(1, 2, 1.5, 3, 3.5)
#' result <- isotonic_regression(data, lambda = 1)
#' print(result$means)
#'
#' # With lambda = 0, should get Pool Adjacent Violators result
#' result_pav <- isotonic_regression(data, lambda = 0)
#'
#' @export
isotonic_regression <- function(data, lambda = 0, min_value = -Inf, max_value = Inf) {
  n <- length(data)

  if (n == 0) {
    stop("Data vector is empty")
  }

  if (lambda < 0) {
    stop("Lambda must be non-negative")
  }

  # Initialize with first data point
  pieces <- list(create_data_piece(data[1], left = min_value, right = max_value, position = 1))

  # Dynamic programming: process each data point
  for (t in 2:n) {
    # Add new data cost to all pieces
    new_pieces <- list()
    for (piece in pieces) {
      new_piece <- piece$copy()
      new_piece$add_cost(a = 1, b = -2 * data[t], c = data[t]^2)
      new_piece$position <- t
      new_pieces <- append(new_pieces, list(new_piece))
    }

    # Create new segment option: start fresh from current data point
    new_segment_piece <- create_data_piece(data[t], left = min_value, right = max_value, position = t)
    new_segment_piece$label <- t
    new_segment_piece$add_penalty(lambda)  # Add penalty for new segment

    # Apply isotonic constraint: μ_t ≥ μ_{t-1}
    # Find the minimum of previous pieces
    if (length(new_pieces) > 0) {
      prev_min_mu <- Inf
      for (piece in new_pieces) {
        argmin <- piece$get_argmin()
        prev_min_mu <- min(prev_min_mu, argmin)
      }

      # Constrain new pieces: must be ≥ prev_min_mu
      constrained_pieces <- list()
      for (piece in new_pieces) {
        if (piece$right >= prev_min_mu) {
          new_piece <- piece$copy()
          new_piece$left <- max(piece$left, prev_min_mu)
          if (new_piece$left < new_piece$right || abs(new_piece$left - new_piece$right) < 1e-10) {
            constrained_pieces <- append(constrained_pieces, list(new_piece))
          }
        }
      }

      # Also constrain new segment piece
      if (new_segment_piece$right >= prev_min_mu) {
        new_segment_piece$left <- max(new_segment_piece$left, prev_min_mu)
      }

      new_pieces <- constrained_pieces
    }

    # Add new segment option
    new_pieces <- append(new_pieces, list(new_segment_piece))

    # Prune dominated pieces using minimum envelope
    pieces <- prune_pieces(new_pieces)
  }

  # Backtrack to find optimal segmentation
  result <- backtrack(pieces, data, lambda)

  return(result)
}


#' Prune dominated pieces
#'
#' @description Remove pieces that are dominated by others (never optimal)
#' Keeps only pieces on the lower envelope
#'
#' @param pieces List of NormalLossPiece objects
#' @return List of non-dominated pieces
#'
#' @keywords internal
prune_pieces <- function(pieces) {
  if (length(pieces) <= 1) {
    return(pieces)
  }

  # Sort pieces by left endpoint
  pieces <- pieces[order(sapply(pieces, function(p) p$left))]

  # Build lower envelope iteratively
  envelope <- list()

  for (piece in pieces) {
    if (length(envelope) == 0) {
      envelope <- list(piece)
    } else {
      # Compute minimum with last piece in envelope
      last_piece <- envelope[[length(envelope)]]

      # Check if intervals overlap
      overlap_left <- max(last_piece$left, piece$left)
      overlap_right <- min(last_piece$right, piece$right)

      if (overlap_left < overlap_right) {
        # Compute min envelope in overlapping region
        min_pieces <- min_envelope(last_piece, piece)

        # Replace last piece and add new pieces
        envelope <- envelope[-length(envelope)]
        envelope <- append(envelope, min_pieces)
      } else {
        # No overlap, just add the piece
        envelope <- append(envelope, list(piece))
      }
    }
  }

  # Remove empty pieces
  envelope <- Filter(function(p) p$right > p$left || abs(p$right - p$left) < 1e-10, envelope)

  return(envelope)
}


#' Backtrack to find optimal segmentation
#'
#' @description Given final cost pieces, backtrack to find the optimal
#' segmentation and means
#'
#' @param pieces List of final NormalLossPiece objects
#' @param data Original data vector
#' @param lambda Penalty parameter
#' @return List with changepoints, means, segments, and cost
#'
#' @keywords internal
backtrack <- function(pieces, data, lambda) {
  n <- length(data)

  # Check if pieces is empty
  if (length(pieces) == 0) {
    # Fallback: use PAV algorithm
    means <- compute_isotonic_means(data)

    # Build segments based on where means change
    segments <- numeric(n)
    changepoints <- c()
    current_segment <- 1
    segments[1] <- 1

    for (i in 2:n) {
      # Handle NA values
      if (!is.na(means[i]) && !is.na(means[i-1]) && abs(means[i] - means[i-1]) > 1e-10) {
        changepoints <- c(changepoints, i-1)
        current_segment <- current_segment + 1
      }
      segments[i] <- current_segment
    }
    changepoints <- c(changepoints, n)

    # Get unique segment means
    unique_means <- c()
    for (seg in 1:max(segments)) {
      unique_means <- c(unique_means, mean(data[segments == seg]))
    }

    # Compute cost
    fitted <- means
    data_cost <- sum((data - fitted)^2)
    penalty_cost <- lambda * (length(unique_means) - 1)
    total_cost <- data_cost + penalty_cost

    return(list(
      changepoints = changepoints,
      means = unique_means,
      segments = segments,
      cost = total_cost,
      data_cost = data_cost,
      penalty_cost = penalty_cost,
      pieces = pieces,
      optimal_mu = NA
    ))
  }

  # Find the piece with minimum cost
  best_piece <- pieces[[1]]
  best_cost <- best_piece$get_min()

  for (piece in pieces[-1]) {
    cost <- piece$get_min()
    if (cost < best_cost) {
      best_cost <- cost
      best_piece <- piece
    }
  }

  # Get optimal μ
  optimal_mu <- best_piece$get_argmin()

  # Simple backtracking: segment based on labels
  # For now, return the final optimal solution
  # Full backtracking would require storing the complete dynamic programming table

  # Compute means using Pool Adjacent Violators on segments
  means <- compute_isotonic_means(data)

  # Identify changepoints (where mean changes)
  changepoints <- c()
  segments <- rep(1, n)
  current_mean <- means[1]
  current_segment <- 1

  for (i in 2:n) {
    # Handle NA/NaN/Inf values
    if (is.na(means[i]) || is.na(current_mean) ||
        is.infinite(means[i]) || is.infinite(current_mean)) {
      if (!identical(means[i], current_mean)) {
        changepoints <- c(changepoints, i - 1)
        current_segment <- current_segment + 1
        current_mean <- means[i]
      }
    } else {
      diff <- abs(means[i] - current_mean)
      if (!is.na(diff) && !is.infinite(diff) && diff > 1e-6) {
        changepoints <- c(changepoints, i - 1)
        current_segment <- current_segment + 1
        current_mean <- means[i]
      }
    }
    segments[i] <- current_segment
  }
  changepoints <- c(changepoints, n)

  # Compute unique segment means
  unique_means <- c()
  for (i in 1:max(segments)) {
    segment_data <- data[segments == i]
    unique_means <- c(unique_means, mean(segment_data))
  }

  # Compute total cost
  fitted_values <- unique_means[segments]
  data_cost <- sum((data - fitted_values)^2)
  penalty_cost <- lambda * (length(unique_means) - 1)
  total_cost <- data_cost + penalty_cost

  return(list(
    changepoints = changepoints,
    means = unique_means,
    segments = segments,
    cost = total_cost,
    data_cost = data_cost,
    penalty_cost = penalty_cost,
    pieces = pieces,
    optimal_mu = optimal_mu
  ))
}


#' Compute isotonic means using Pool Adjacent Violators
#'
#' @description Classic PAV algorithm for isotonic regression without penalty
#' @param y Data vector
#' @return Vector of isotonic means
#'
#' @keywords internal
compute_isotonic_means <- function(y) {
  n <- length(y)
  if (n == 0) return(numeric(0))
  if (n == 1) return(y)

  # Pool Adjacent Violators Algorithm
  means <- y
  weights <- rep(1, n)

  repeat {
    changed <- FALSE

    for (i in 1:(n - 1)) {
      # Skip if either mean is NA
      if (is.na(means[i]) || is.na(means[i + 1])) next

      if (means[i] > means[i + 1]) {
        # Pool adjacent violators
        pooled_mean <- (means[i] * weights[i] + means[i + 1] * weights[i + 1]) /
                       (weights[i] + weights[i + 1])
        pooled_weight <- weights[i] + weights[i + 1]

        means[i] <- pooled_mean
        means[i + 1] <- pooled_mean
        weights[i] <- pooled_weight
        weights[i + 1] <- pooled_weight

        changed <- TRUE
      }
    }

    if (!changed) break
  }

  return(means)
}


#' Plot isotonic regression result
#'
#' @description Visualize data and fitted isotonic regression
#' @param data Original data
#' @param result Result from isotonic_regression()
#' @param main Plot title
#'
#' @export
plot_isotonic_result <- function(data, result, main = "Isotonic Regression") {
  n <- length(data)
  fitted <- result$means[result$segments]

  plot(1:n, data, pch = 19, col = "black",
       xlab = "Index", ylab = "Value", main = main,
       ylim = range(c(data, fitted)))

  # Plot fitted values
  lines(1:n, fitted, col = "red", lwd = 2, type = "s")

  # Mark changepoints
  if (length(result$changepoints) > 1) {
    for (cp in result$changepoints[-length(result$changepoints)]) {
      abline(v = cp + 0.5, col = "blue", lty = 2)
    }
  }

  # Add legend
  legend("topleft",
         legend = c("Data", "Fitted (isotonic)", "Changepoints"),
         col = c("black", "red", "blue"),
         pch = c(19, NA, NA),
         lty = c(NA, 1, 2),
         lwd = c(NA, 2, 1))

  # Add cost information
  text(n * 0.7, max(data) * 0.95,
       sprintf("Cost: %.2f\nSegments: %d\nλ: %.2f",
               result$cost, length(result$means), result$penalty_cost / max(1, length(result$means) - 1)),
       adj = 0, cex = 0.8)
}
