#' Compute the optimal imputation minimizing a generalized dissimilarity
#'
#' Solves the quadratic program
#' \deqn{
#' \min_{x \in \mathbb{R}^n}
#'   \sum_{S \subseteq N} \gamma_S \left(v(S) - \sum_{i \in S} x_i\right)^2
#'   + \sum_{S \subseteq N} \beta_S \left(v(S) - \sum_{i \in S} x_i\right)
#' }
#' Note that the solution for \eqn{\gamma_S} equal to the Shapley weights and \eqn{\beta_S = 0}
#' recovers the Shapley value.
#'
#' @param v Numeric vector of length \eqn{(2^n-1)}, the game values (\eqn{\gamma_S}) for all non‑empty coalitions.
#' @param gamma Numeric weight vector of the same length as \code{v}, (\eqn{\gamma_S}) in the objective.
#' @param beta  Numeric bias vector of the same length as \code{v}, (\eqn{\beta_S}) in the objective. Defaults to zero.
#'
#' @return Numeric vector \eqn{x^* \in \mathbb{R}^n}, the optimal imputation minimizing the given objective.
#'
#' @examples
#' v      <- c(0, 0.8, 0.8, 0, 0.8, 0, 1)
#' n      <- log2(length(v) + 1)
#' gamma  <- get_shapley_gamma(n, last_val = v[length(v)])
#' beta   <- rep(0, length(v))
#' psi    <- get_psi_allocation(v, gamma = gamma, beta = beta)
#' shap   <- CoopGame::shapleyValue(v)
#'
#' @importFrom quadprog solve.QP
#' @importFrom CoopGame createBitMatrix shapleyValue
#' @export
get_psi_allocation <- function(v,
                               gamma = rep(1, length(v)),
                               beta  = rep(0, length(v))) {
  # Validate dimensions
  n <- as.integer(log2(length(v) + 1))
  if (!isTRUE(all.equal(2^n - 1, length(v)))) {
    stop("Length of v must be 2^n - 1 for some integer n.")
  }
  if (length(gamma) != length(v) || length(beta) != length(v)) {
    stop("gamma and beta must have the same length as v.")
  }

  # Precompute coalition membership matrix (2^n-1 x n)
  M <- createBitMatrix(n)[, 1:n]

  # Build Hessian Dmat (n x n) and linear term dvec
  # Dmat[i,j] = 2 * sum_{S: i,j in S} gamma_S
  # dvec[i]   = sum_{S: i in S} (2*gamma_S*v(S) + beta_S)
  Dmat <- matrix(0, nrow = n, ncol = n)
  dvec <- numeric(n)

  # Compute contributions by matrix multiplication
  # For off-diagonals: M^T %*% diag(gamma) %*% M gives counts of gamma_S for each pair (i,j)
  G  <- diag(gamma)
  mat <- t(M) %*% G %*% M
  Dmat[] <- 2 * mat      # multiply entire matrix by 2

  # Diagonal got same formula from mat[i,i]
  # Build dvec: t(M) %*% (2*gamma*v + beta)
  delta <- 2 * gamma * v + beta
  dvec  <- as.numeric(t(M) %*% delta)

  # Equality constraint: sum(x) == v(N)
  Amat  <- matrix(1, nrow = n, ncol = 1)
  bvec  <- v[length(v)]

  # Solve QP: minimize 1/2 x^T D x - d^T x subject to A^T x = b
  # To enforce equality, set meq = 1

  result <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1L)

  return(result$solution)
}
