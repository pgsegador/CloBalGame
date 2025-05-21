#' Finds the closest balanced game (quadratic formulation)
#'
#' This implements the Goldfarb–Idnani dual algorithm (via **quadprog**) to solve
#' the following constrained weighted least-squares problem:
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
#' @param v A numeric vector of length \eqn{2^n - 1} representing a transferable utility (TU) game.
#' @param gamma A non-negative weight vector of the same length as \code{v}. Defaults to uniform weights associated to Euclidean distance.
#' @param vN Optional numeric value to fix \eqn{v^*(N)} in the optimization. By default, this value is set to match the value of the original game \code{vN=v[length(v)]}. A valid value is \code{vN = NULL}; in this case, the value of the balanced game for the grand coalition, \eqn{v^*(N)}, is not fixed and is instead optimized within the optimization problem.
#' @param positive_game Logical; if \code{TRUE}, enforces \eqn{v^*(S) \ge 0} for all coalitions.
#' @param epsilon A small positive scalar (default \code{1e-8}) added to the diagonal of the Hessian to ensure positive definiteness and numerical stability.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{v_star}}{Numeric vector representing the closest balanced game \eqn{v^*}.}
#'   \item{\code{x_star}}{Numeric vector representing the associated core imputation \eqn{x^*}.}
#' }
#'
#' @examples
#' v <- c(4,2,1,5,6,3,4)
#' closest_balanced_game_quad(v)
#'
#' @importFrom quadprog solve.QP
#' @importFrom CoopGame createBitMatrix
#' @export
closest_balanced_game_quad <- function(
    v,
    gamma = rep(1, length(v)),
    vN = v[length(v)],
    positive_game = FALSE,
    epsilon = 1e-8
) {
  # Number of players
  n <- as.integer(log2(length(v) + 1))
  if (!isTRUE(all.equal(2^n - 1, length(v))))
    stop("Length of v must be 2^n - 1 for some integer n.")
  if (length(gamma) != length(v))
    stop("gamma must have the same length as v.")
  if (positive_game && !is.null(vN) && vN < 0)
    stop("For positive_game = TRUE and vN != NULL, vN must be positive (vN >= 0).")


  # 0) Build quadratic term Dmat and linear term dvec
  Dmat <- diag(c(2 * gamma, rep(epsilon, n)))
  dvec <- c(2 * gamma * v, rep(0, n))

  ## Matrix defining the constraints

  # Precompute sizes
  num_coalitions <- length(v)        # 2^n - 1
  num_players    <- n                # as derived earlier
  total_vars     <- num_coalitions + num_players  # for Dmat padding


  # 1) Build the coalition‐membership part
  Mat_S <- createBitMatrix(n)[, 1:num_players]

  # 2) Combine with the “–I” part for coalition‐rationality:  x(S) ≥ v(S)
  Mat_S <- cbind(-diag(num_coalitions), Mat_S)

  # 3) Move the grand‐coalition row to the top so equalities come first
  Mat_S <- rbind(Mat_S[num_coalitions,],Mat_S)
  Mat_S <- Mat_S[-(num_coalitions+1),]

  # Default: only one equality (efficiency)
  meq <- 1L

  # 4) If vN is specified, add its exact‐value constraint
  if (!is.null(vN)) {
    eq_row     <- integer(ncol(Mat_S))
    eq_row[num_coalitions] <- 1         # picks out x(N)
    Mat_S      <- rbind(eq_row, Mat_S)
    meq        <- 2L                     # now two equalities at top
  }

  # 5) If positivity is required, append x ≥ 0 constraints
  if (positive_game) {
    pos_block <- cbind(diag(num_coalitions), matrix(0, num_coalitions, num_players))
    Mat_S     <- rbind(Mat_S, pos_block)
  }

  # 6) Build the right‐hand side vector
  bvec <- rep(0, nrow(Mat_S))
  if (!is.null(vN)) {
    bvec[1] <- vN  # fix x(N) to vN
  }

  # 7) Transpose for quadprog:  solve.QP( D, d, G = t(Mat_S), bvec, meq )
  Mat_S <- t(Mat_S)

  qp <- solve.QP(Dmat, dvec, Mat_S, bvec, meq = meq)

  # 8) Return v_star, x_star and opt_value

  # Extract solutions
  v_star <- qp$solution[seq_len(num_coalitions)]
  x_star <- qp$solution[(num_coalitions+1):total_vars]

  out <- list(
    v_star = v_star,
    x_star = x_star
  )

  return(out)

}
