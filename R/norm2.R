#' Norm2
#'
#' @param x
#'
#' @return Returns the l2-norm of a vector
#' @export
#'
#' @examples
#' norm1 (rnorm(10))


norm2 <- function(x) {
  sqrt(sum(x ^ 2))
}
