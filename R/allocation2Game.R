#' Converts an allocation to its associated additive game
#'
#' Given an allocation vector \eqn{x \in \mathbb{R}^n}, this function constructs
#' the additive game \eqn{v} such that \eqn{v(S) = \sum_{i \in S} x_i} for all coalitions \eqn{S}.
#'
#' @param x A numeric vector representing an allocation (payoff for each player).
#'
#' @return A numeric vector of length \eqn{2^n - 1} representing the additive game associated with \code{x}.
#'
#' @examples
#' x <- c(0.1, 0.4, 0.5)
#' v <- allocation2Game(x)
#' print(v)
#'
#' @importFrom CoopGame createBitMatrix
#'
#' @export
allocation2Game <- function(x) {

  n <- length(x)
  # Create the binary matrix for all coalitions and compute v(S) = sum_{i in S} x_i
  v <- as.numeric(createBitMatrix(n)[, 1:n] %*% x)

  return(v)
}
