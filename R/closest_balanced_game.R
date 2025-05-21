#' Compute the Closest Balanced Game via Quadratic or Iterative Methods
#'
#' Calculates the closest balanced game \eqn{v^*} and associated core imputation \eqn{x^*} for a cooperative TU game.
#' In other words, it solves the following problem:
#' \deqn{
#' \min_{v^*,\, x^*} \sum_{S} \gamma_S\, (v(S) - v^*(S))^2
#' }
#' subject to the constraints:
#' \deqn{
#' \sum_{i \in N} x^*_i = v^*(N) \quad \text{(efficiency)}
#' }
#' \deqn{
#' \sum_{i \in S} x^*_i \ge v^*(S) \quad \text{for all } S \subsetneq N \quad \text{(coalitional rationality)}
#' }
#' Optionally, the following constraints can also be imposed:
#' \itemize{
#'   \item \eqn{v^*(S) \ge 0} for all coalitions (\code{positive_game}).
#'   \item \eqn{v^*(N) = v_N} to fix the value of the grand coalition (\code{vN}). By default, this value is set to match the value of the original game \eqn{v_N = v(N)}.
#' }
#'
#'
#' @param v Numeric vector of length \eqn{2^n - 1}, game values for each non-empty coalition.
#' @param gamma Numeric vector of same length as \code{v}, weights for each coalition.
#' @param vN Optional numeric value to fix \eqn{v^*(N)} in the optimization. By default, this is set to the grand coalition value \code{v[length(v)]}. A valid value is \code{vN = NULL}, in which case \eqn{v^*(N)} is optimized as part of the problem.
#' @param positive_game Logical; if \code{TRUE}, enforces non-negativity on the solution (all \eqn{v^*(S) \ge 0}).
#' @param method Character; either \code{'iter'} (iterative) or \code{'quad'} (quadratic). Default is \code{'iter'}.
#' @param tol Numeric tolerance for convergence in the iterative method (default = 1e-22).
#' @param max_iter_seed Integer, maximum iterations before reseeding in the iterative method (default = 2*n).
#' @param prob_seed Numeric in (0,1), probability to include each coalition when reseeding (default = 0.8). Only used if \code{method = 'iter'}.
#' @param epsilon A small positive scalar (default = 1e-8) added to the diagonal of the Hessian for numerical stability. Only used if \code{method = 'quad'}.
#'
#' @return A list with components:
#' \describe{
#'   \item{v_star}{Numeric vector of length \eqn{2^n - 1}, the closest balanced game values.}
#'   \item{x_star}{Numeric vector of length \eqn{n}, the associated imputation.}
#' }
#'
#' @details
#' Depending on \code{method}, this function either:
#' \itemize{
#'   \item \code{'iter'}: uses an active-set iterative algorithm with probabilistic reseeding. The \code{"iter"} (iterative) method, if it terminates, guarantees the optimal value of the closest balanced game. However, it may get stuck in a loop, especially for small values of \eqn{n} (e.g., \eqn{n = 3} or \eqn{n = 4}). It is generally faster and more efficient for larger games (\eqn{n > 10}).
#'   \item \code{'quad'}: formulates and solves a quadratic program enforcing efficiency and optional positivity. The \code{"quad"} (quadratic programming) method does not suffer from infinite loops and always provides a result. However, it is not recommended for large games (\eqn{n > 10}) due to poor computational efficiency.
#' }
#'
#' @examples
#' v <- c(4,2,1,5,6,3,4)
#' # Iterative method (default)
#' closest_balanced_game(v, method = 'iter')
#' # Quadratic method with positivity
#' closest_balanced_game(v, method = 'quad', positive_game = TRUE)
#'
#' @importFrom MASS ginv
#' @importFrom quadprog solve.QP
#' @importFrom CoopGame createBitMatrix
#' @importFrom stats runif
#' @export
closest_balanced_game <- function(
    v,
    gamma = rep(1, length(v)),
    vN = v[length(v)],
    positive_game = FALSE,
    method = 'iter',
    tol = 1e-22,
    max_iter_seed = 2 * as.integer(log2(length(v) + 1)),
    prob_seed = 0.8,
    epsilon = 1e-8
) {

  n <- as.integer(log2(length(v) + 1))
  if (method == "quad" && n > 10) {
    warning("Method 'quad' is not recommended for large games (n > 10). Consider using method = 'iter' instead.")
  }

  if (method == 'iter') {
    res <- closest_balanced_game_iter(
      v = v,
      gamma = gamma,
      vN = vN,
      positive_game = positive_game,
      tol = tol,
      max_iter_seed = max_iter_seed,
      prob_seed = prob_seed,
      epsilon = epsilon
    )
    return(res)

  } else if (method == 'quad') {
    res_quad <- closest_balanced_game_quad(
      v = v,
      gamma = gamma,
      vN = vN,
      positive_game = positive_game,
      epsilon = epsilon
    )
    return(res_quad)
  } else {
    stop("Unknown method. Choose 'iter' or 'quad'.")
  }
}

