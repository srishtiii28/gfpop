# NormalLossPiece Class for Regularized Isotonic Regression
# Implements piecewise quadratic cost functions for normal (Gaussian) loss

#' NormalLossPiece R6 Class
#'
#' @description Represents a piece of a piecewise quadratic cost function.
#' Each piece is valid over an interval [left, right] and has the form:
#' cost(μ) = A*μ² + B*μ + C
#'
#' This class is essential for implementing regularized isotonic regression
#' using dynamic programming with piecewise cost functions.
#'
#' @field A Numeric. Quadratic coefficient (always positive for normal loss)
#' @field B Numeric. Linear coefficient
#' @field C Numeric. Constant term
#' @field left Numeric. Left endpoint of the interval where this piece is valid
#' @field right Numeric. Right endpoint of the interval where this piece is valid
#' @field label Integer. Segment label for backtracking
#' @field position Integer. Data position where this piece originates
#'
#' @examples
#' # Create a piece for data point y = 5
#' # Cost = (y - μ)² = μ² - 2*5*μ + 25
#' piece <- NormalLossPiece$new(A = 1, B = -10, C = 25, left = -Inf, right = Inf)
#' piece$eval(5)  # Should be 0 (exact fit)
#' piece$get_argmin()  # Should be 5
#'
#' @export
library(R6)

NormalLossPiece <- R6Class("NormalLossPiece",
  public = list(
    A = NULL,
    B = NULL,
    C = NULL,
    left = NULL,
    right = NULL,
    label = NULL,
    position = NULL,

    #' @description Initialize a new NormalLossPiece
    #' @param A Quadratic coefficient
    #' @param B Linear coefficient
    #' @param C Constant term
    #' @param left Left endpoint of interval (default: -Inf)
    #' @param right Right endpoint of interval (default: Inf)
    #' @param label Segment label (default: 0)
    #' @param position Data position (default: 0)
    initialize = function(A = 0, B = 0, C = 0, left = -Inf, right = Inf, label = 0, position = 0) {
      self$A <- A
      self$B <- B
      self$C <- C
      self$left <- left
      self$right <- right
      self$label <- label
      self$position <- position
    },

    #' @description Evaluate cost function at a given μ value
    #' @param mu The parameter value to evaluate
    #' @return The cost value
    eval = function(mu) {
      if (mu < self$left || mu > self$right) {
        return(Inf)
      }
      return(self$A * mu^2 + self$B * mu + self$C)
    },

    #' @description Get the argmin of the quadratic over its valid interval
    #' @return The μ value that minimizes the cost within [left, right]
    get_argmin = function() {
      if (abs(self$A) < 1e-10) {
        # Nearly linear, return boundary
        if (abs(self$B) < 1e-10) {
          return((self$left + self$right) / 2)
        }
        return(if (self$B < 0) self$right else self$left)
      }

      # Quadratic case: argmin = -B / (2*A)
      unconstrained_argmin <- -self$B / (2 * self$A)

      # Handle NA/NaN/Inf cases
      if (is.na(unconstrained_argmin) || is.infinite(unconstrained_argmin)) {
        # If argmin is invalid, return midpoint of interval
        if (is.infinite(self$left) && is.infinite(self$right)) {
          return(0)
        } else if (is.infinite(self$left)) {
          return(self$right)
        } else if (is.infinite(self$right)) {
          return(self$left)
        } else {
          return((self$left + self$right) / 2)
        }
      }

      # Clip to interval
      if (unconstrained_argmin < self$left) {
        return(self$left)
      } else if (unconstrained_argmin > self$right) {
        return(self$right)
      } else {
        return(unconstrained_argmin)
      }
    },

    #' @description Get the minimum cost value over the interval
    #' @return The minimum cost
    get_min = function() {
      argmin <- self$get_argmin()
      return(self$eval(argmin))
    },

    #' @description Add another quadratic cost to this piece
    #' @param a Quadratic coefficient to add
    #' @param b Linear coefficient to add
    #' @param c Constant term to add
    add_cost = function(a, b, c) {
      self$A <- self$A + a
      self$B <- self$B + b
      self$C <- self$C + c
    },

    #' @description Add a constant penalty
    #' @param penalty Penalty value to add
    add_penalty = function(penalty) {
      self$C <- self$C + penalty
    },

    #' @description Create a copy of this piece
    #' @return A new NormalLossPiece object with the same values
    copy = function() {
      NormalLossPiece$new(
        A = self$A, B = self$B, C = self$C,
        left = self$left, right = self$right,
        label = self$label, position = self$position
      )
    },

    #' @description Get cost from data point (y - μ)²
    #' @param y Data value
    #' @return List with A, B, C coefficients
    from_data = function(y) {
      # (y - μ)² = μ² - 2*y*μ + y²
      list(A = 1, B = -2 * y, C = y^2)
    },

    #' @description Print piece information
    print = function() {
      cat(sprintf("NormalLossPiece: [%.2f, %.2f]\n", self$left, self$right))
      cat(sprintf("  Cost: %.4f*μ² + %.4f*μ + %.4f\n", self$A, self$B, self$C))
      cat(sprintf("  Argmin: %.4f, Min: %.4f\n", self$get_argmin(), self$get_min()))
      cat(sprintf("  Label: %d, Position: %d\n", self$label, self$position))
    }
  )
)


#' Create NormalLossPiece from data point
#'
#' @description Helper function to create a NormalLossPiece for a single data point
#' with loss (y - μ)²
#'
#' @param y Data value
#' @param left Left endpoint of interval
#' @param right Right endpoint of interval
#' @param position Data position
#' @return A new NormalLossPiece object
#'
#' @examples
#' piece <- create_data_piece(y = 5, position = 1)
#' piece$get_argmin()  # Should be 5
#'
#' @export
create_data_piece <- function(y, left = -Inf, right = Inf, position = 0) {
  # Cost = (y - μ)² = μ² - 2*y*μ + y²
  NormalLossPiece$new(
    A = 1,
    B = -2 * y,
    C = y^2,
    left = left,
    right = right,
    label = position,
    position = position
  )
}


#' Compute minimum envelope of two pieces
#'
#' @description Given two NormalLossPiece objects, compute their minimum envelope.
#' This finds the intervals where each piece has lower cost and returns
#' a list of pieces representing min(piece1, piece2).
#'
#' @param piece1 First NormalLossPiece
#' @param piece2 Second NormalLossPiece
#' @return List of NormalLossPiece objects representing the minimum envelope
#'
#' @export
min_envelope <- function(piece1, piece2) {
  # Find where piece1 == piece2
  # A1*μ² + B1*μ + C1 = A2*μ² + B2*μ + C2
  # (A1-A2)*μ² + (B1-B2)*μ + (C1-C2) = 0

  dA <- piece1$A - piece2$A
  dB <- piece1$B - piece2$B
  dC <- piece1$C - piece2$C

  # Find intersection points
  roots <- c()
  if (abs(dA) < 1e-10) {
    # Linear case
    if (abs(dB) > 1e-10) {
      roots <- c(-dC / dB)
    }
  } else {
    # Quadratic case
    discriminant <- dB^2 - 4 * dA * dC
    if (!is.na(discriminant) && !is.infinite(discriminant) && discriminant >= 0) {
      sqrt_disc <- sqrt(discriminant)
      roots <- c(
        (-dB - sqrt_disc) / (2 * dA),
        (-dB + sqrt_disc) / (2 * dA)
      )
    }
  }

  # Filter roots to valid interval
  interval_left <- max(piece1$left, piece2$left)
  interval_right <- min(piece1$right, piece2$right)

  # Check if intervals overlap
  if (interval_left >= interval_right) {
    return(list())
  }

  roots <- roots[!is.na(roots) & !is.infinite(roots)]
  roots <- roots[roots > interval_left & roots < interval_right]
  roots <- sort(unique(roots))

  # Build envelope pieces
  result <- list()
  breakpoints <- c(interval_left, roots, interval_right)

  for (i in 1:(length(breakpoints) - 1)) {
    left <- breakpoints[i]
    right <- breakpoints[i + 1]

    if (is.infinite(left) || is.infinite(right)) next
    if (right - left < 1e-10) next

    # Test which piece is lower at midpoint
    mid <- (left + right) / 2
    eval1 <- piece1$eval(mid)
    eval2 <- piece2$eval(mid)

    # Handle NA/NaN cases
    if (is.na(eval1) && is.na(eval2)) next
    if (is.na(eval1)) {
      new_piece <- piece2$copy()
    } else if (is.na(eval2)) {
      new_piece <- piece1$copy()
    } else if (eval1 <= eval2) {
      new_piece <- piece1$copy()
    } else {
      new_piece <- piece2$copy()
    }

    new_piece$left <- left
    new_piece$right <- right
    result <- append(result, list(new_piece))
  }

  return(result)
}


#' Visualize cost function pieces
#'
#' @description Plot piecewise quadratic cost function
#' @param pieces List of NormalLossPiece objects
#' @param title Plot title
#' @param xlim X-axis limits
#' @param add_points Add points at breakpoints
#'
#' @export
plot_pieces <- function(pieces, title = "Piecewise Quadratic Cost", xlim = NULL, add_points = TRUE) {
  if (length(pieces) == 0) {
    stop("No pieces to plot")
  }

  # Determine plot range
  if (is.null(xlim)) {
    all_lefts <- sapply(pieces, function(p) p$left)
    all_rights <- sapply(pieces, function(p) p$right)
    all_lefts[is.infinite(all_lefts)] <- NA
    all_rights[is.infinite(all_rights)] <- NA

    xlim <- c(
      min(all_lefts, na.rm = TRUE) - 1,
      max(all_rights, na.rm = TRUE) + 1
    )
  }

  # Generate plot points
  x_vals <- seq(xlim[1], xlim[2], length.out = 500)
  y_vals <- numeric(length(x_vals))

  for (i in seq_along(x_vals)) {
    x <- x_vals[i]
    # Find which piece contains x
    y <- Inf
    for (piece in pieces) {
      if (x >= piece$left && x <= piece$right) {
        y <- min(y, piece$eval(x))
      }
    }
    y_vals[i] <- y
  }

  # Remove infinite values for plotting
  valid <- is.finite(y_vals)
  x_vals <- x_vals[valid]
  y_vals <- y_vals[valid]

  plot(x_vals, y_vals, type = "l", lwd = 2, col = "blue",
       xlab = expression(mu), ylab = "Cost",
       main = title)

  if (add_points) {
    # Add breakpoints
    breakpoints <- c()
    for (piece in pieces) {
      if (is.finite(piece$left)) breakpoints <- c(breakpoints, piece$left)
      if (is.finite(piece$right)) breakpoints <- c(breakpoints, piece$right)
    }
    breakpoints <- unique(breakpoints)

    for (bp in breakpoints) {
      if (bp >= xlim[1] && bp <= xlim[2]) {
        y <- Inf
        for (piece in pieces) {
          if (bp >= piece$left && bp <= piece$right) {
            y <- min(y, piece$eval(bp))
          }
        }
        if (is.finite(y)) {
          points(bp, y, pch = 19, col = "red", cex = 1.2)
        }
      }
    }
  }

  grid()
}
