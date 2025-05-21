#' Generate Coalition Labels for a Cooperative Game
#'
#' This function generates character labels representing all possible non-empty coalitions
#' of players in a cooperative game with `n` players. Each coalition is represented as a
#' concatenated string of player indices (e.g., "12", "13", "123").
#'
#' @param n Integer. The number of players in the game.
#'
#' @return A character vector with all non-empty coalitions represented as strings.
#' The coalitions are ordered first by size and then lexicographically.
#'
#' @examples
#' get_index_set(3)
#' # Returns: "1" "2" "3" "12" "13" "23" "123"
#'
#' @importFrom utils combn
#'
#' @export
get_index_set <- function(n) {
  # Initialize character vector to store coalition labels
  coalitions <- character()

  # Generate all combinations of players for coalition sizes 1 to n
  for (i in 1:n) {
    combinations <- combn(1:n, i, simplify = FALSE)
    labels <- vapply(combinations, function(x) paste(x, collapse = ""), character(1))
    coalitions <- c(coalitions, labels)
  }

  return(coalitions)
}
