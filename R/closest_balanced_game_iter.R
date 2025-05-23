#' Iteratively Compute the Closest Balanced Imputation with Probabilistic Reseeding
#'
#' Given a cooperative game vector \code{v} and weights \code{gamma}, this function
#' iteratively finds a core imputation of \eqn{v^*}, and uses this imputation to deduce the value of the closest balanced game \eqn{v^*}.
#'
#' @param v Numeric vector of length \eqn{2^n - 1}, game values for each non-empty coalition.
#' @param gamma Numeric vector of same length as \code{v}, weights for each coalition.
#' @param vN Optional numeric value to fix \eqn{v^*(N)} in the optimization. By default, this value is set to match the value of the original game \code{vN=v[length(v)]}. A valid value is \code{vN = NULL}; in this case, the value of the balanced game for the grand coalition, \eqn{v^*(N)}, is not fixed and is instead optimized within the optimization problem.
#' @param positive_game Logical; if \code{TRUE}, enforces \eqn{v^*(S) \ge 0} for all coalitions.
#' @param tol Numeric tolerance for feasibility comparison in each iteration (default = 1e-22).
#' @param max_iter_seed Integer, maximum iterations before reseeding (default = 2*n).
#' @param prob_seed Numeric in (0,1), probability to include each coalition when reseeding (default = 0.8).
#' @param epsilon A small positive scalar (default \code{1e-8}) added to the diagonal of the Hessian to ensure positive definiteness and numerical stability.
#'
#' @return Numeric vector \code{x_star} of length \code{n}, the final imputation.
#'
#' @details
#' - At most \code{max_iter_seed} iterations are allowed before a random reseed of the index set.
#' - On reseed, a Bernoulli trial with probability \code{prob_seed} determines which coalitions to include.
#' - This helps avoid infinite loops when convergence stalls.
#'
#' @examples
#' v <- c(4,2,1,5,6,3,4)
#' closest_balanced_game_iter(v)
#'
#' @import MASS
#' @importFrom quadprog solve.QP
#' @importFrom CoopGame createBitMatrix
#' @importFrom stats runif
#'
#' @export
closest_balanced_game_iter <- function(
    v,
    gamma = rep(1, length(v)),
    vN = v[length(v)],
    positive_game = FALSE,
    tol = 1e-22,
    max_iter_seed = 2 * as.integer(log2(length(v) + 1)),
    prob_seed = 0.8,
    epsilon = 1e-8
) {

  # Determine number of players and coalitions
  n <- as.integer(log2(length(v) + 1))
  m <- length(v)

  if (2^n - 1 != m)
    stop("Length of v must be 2^n - 1 for some integer n.")
  if (length(gamma) != m)
    stop("gamma must have the same length as v.")
  if (positive_game && !is.null(vN) && vN < 0)
    stop("For positive_game = TRUE and vN != NULL, vN must be positive (vN >= 0).")

  # Binary coalition membership matrix (m x n)
  B <- CoopGame::createBitMatrix(n)[, 1:n, drop = FALSE]

  # Initialize active set, iteration counter, and reseed counter
  A_prev <- NULL
  iter <- 0L

  repeat {
    iter <- iter + 1L

    # If reached max_iter_seed without convergence, reseed active set randomly
    if (iter > max_iter_seed) {
      A_prev <- NULL
      bern <- runif(m) < prob_seed
      idx <- which(bern)
      iter <- 1L
    } else {
      # Choose indices: all on first pass, else previous active set
      idx <- if (is.null(A_prev)) seq_len(m) else A_prev
    }

    # Subset matrices and vectors
    B_sub <- B[idx, , drop = FALSE]
    g_sub <- gamma[idx]
    v_sub <- v[idx]

    ## CASE 1
    if (!is.null(vN) && !positive_game){

      # Weighted moment matrix and RHS
      G_B <- sweep(B_sub, 1, g_sub, "*")        # each row * gamma
      M_mat <- t(B_sub) %*% G_B                   # n x n
      N_vec <- as.numeric(t(B_sub) %*% (g_sub * v_sub))


      # Build equations: differences for players 1..(n-1)
      if (n > 1) {
        A_top <- M_mat[1:(n - 1), , drop = FALSE] - M_mat[2:n, , drop = FALSE]
        b_top <- N_vec[1:(n - 1)] - N_vec[2:n]
      } else {
        A_top <- matrix(0, 0, n)
        b_top <- numeric(0)
      }

      # Efficiency equation: sum(x) = vN
      A_last <- rep(1, n)
      b_last <- vN

      # Combine to full system
      A_mat <- rbind(A_top, A_last)
      b_vec <- c(b_top, b_last)

      # Solve linear system by pseudoinverse
      x_star <- as.vector(MASS::ginv(A_mat) %*% b_vec)

      ## CASE 2
    } else if (is.null(vN) && !positive_game){

      # Weighted moment matrix and RHS
      G_B <- sweep(B_sub, 1, g_sub, "*")        # each row * gamma
      M_mat <- t(B_sub) %*% G_B                   # n x n
      N_vec <- as.numeric(t(B_sub) %*% (g_sub * v_sub))

      # Solve the system M_mat %*% x = N_vec
      x_star <- as.vector(MASS::ginv(M_mat) %*% N_vec)

      ## CASE 3
    } else if (!is.null(vN) && positive_game){

      # 0) Quadratic terms
      Dmat <- 2 * t(B_sub) %*% diag(g_sub ,nrow = length(g_sub),ncol = length(g_sub)) %*% B_sub + diag(epsilon, n)  # regularización
      dvec <- 2 * t(B_sub) %*% (g_sub * v_sub)  # linear

      # 1) Efficiency constraint: x(N) = v(N)
      A_eq <- matrix(1, nrow = 1, ncol = n)
      b_eq <- vN  # value of the grand coalition

      # 2) Non-negativity constraints: x_i ≥ 0
      A_ineq <- diag(n)
      b_ineq <- rep(0, n)

      # 3) Combine constraints
      Amat <- t(rbind(A_eq, A_ineq))  # in solve.QP, Amat should be transposed
      bvec <- c(b_eq, b_ineq)
      meq <- 1

      # 4) solve
      x_star <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)$solution

      # CASE 4
    } else if (is.null(vN) && positive_game){

      # 0) Quadratic terms
      Dmat <- 2 * t(B_sub) %*% diag(g_sub ,nrow = length(g_sub),ncol = length(g_sub)) %*% B_sub + diag(epsilon, n)  # regularización
      dvec <- 2 * t(B_sub) %*% (g_sub * v_sub)  # linear

      # 1) Non-negativity constraints: x_i ≥ 0
      A_ineq <- diag(n)
      b_ineq <- rep(0, n)

      # 2) Prepare input
      Amat <- t(A_ineq)  # transposed
      bvec <- b_ineq
      meq <- 0  # no equality constraints

      # 3) solve
      x_star <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)$solution
    }

    # Compute induced game from x_star
    v_alloc <- as.numeric(B %*% x_star)

    # Update active set

    if (!is.null(vN)){

      A_current <- setdiff(which(v >= (v_alloc - tol)), m)

    }else{

      A_current <- c(which(v >= (v_alloc - tol)), m)
    }

    # Check convergence: active set stable
    if (!is.null(A_prev) && identical(A_current, A_prev)) {
      break
    }

    # Prepare for next iteration
    A_prev <- A_current
  }

  # Output

  v_star = pmin(v_alloc,v)

  if (positive_game){ v_star <- pmax(v_star,0) }

  if (is.null(vN)){

    v_star[m] <- v_alloc[m]

  } else {

    v_star[m] <- vN

  }


  out <- list(
    v_star = v_star,
    x_star = x_star
  )

  return(out)

}

