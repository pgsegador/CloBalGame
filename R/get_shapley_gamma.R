#' Compute Gamma Weights for the Shapley-Based Metric
#'
#' This function returns the vector of weights \eqn{\gamma_S} associated with a special
#' metric used in cooperative game theory. This metric is minimized by the Shapley value
#' and is defined as: \deqn{\sqrt{\sum_S \gamma_S (v(S) - x(S))^2}}
#' where the sum runs over all coalitions \eqn{S}.
#'
#' @param n Integer. The number of players in the cooperative game.
#' @param last_val Numeric (default = 1). A custom value for the gamma weight of the grand coalition.
#'
#' @return A numeric vector of length \eqn{2^n - 1} with weights \eqn{\gamma_S} for each coalition.
#'
#' @details
#' These weights are derived from the characterization of the Shapley value as the unique
#' allocation that minimizes a weighted Euclidean distance from the game vector.
#'
#' The value of \eqn{\gamma_S} depends on the cardinality \eqn{|S|} of each coalition, and is given by:
#'
#' \deqn{\gamma_S = \frac{1}{\binom{n-2}{|S|-1}}}
#'
#' The weight for the grand coalition (the last entry) is set to `last_val`.
#'
#' @examples
#' get_shapley_gamma(3)
#'
#' @export
get_shapley_gamma <- function(n, last_val = 1) {
  # Compute number of coalitions of each size (1 to n)
  coalition_sizes <- choose(n, 1:n)

  # Repeat each size according to number of coalitions of that size
  cardinalities <- rep(seq_along(coalition_sizes), times = coalition_sizes)

  # Compute gamma weights
  gamma <- 1 / choose(n - 2, cardinalities - 1)

  # Assign user-defined value to the grand coalition
  gamma[length(gamma)] <- last_val

  return(gamma)
}
