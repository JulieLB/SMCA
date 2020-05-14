#' Create a vector of B12 inter B2 projection
#'
#' @param X A numeric data matrix issued from a qualitative data
#' @param G A partition of the matrix columns
#'
#' @return Returns the vector v containing the norm2 of each group.
#' @export
#'
#' @examples
#'G <- list(1, c(2, 3), c(4, 5))
#'a <- rnorm(5)
#'v_subvector(a, G)
#'

foo_tmp <- function(g, x) norm2(x[g])

v_subvector <- function(x, G) {

  v <- sapply(G, foo_tmp, x = x)
  # v <- c()
  #
  # if (is.vector(x)) {
  #   for (i in 1:length(G)) {
  #     iota <- G[[i]]
  #     v <- c(v, norm2(x[iota]))
  #   }
  # }

  return(v)
}
