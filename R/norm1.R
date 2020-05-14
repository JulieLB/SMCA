#'Norm1
#' @description  This function returns the l1 norm of a vector.
#' @param vect a vector
#'
#' @return  the l1-norm of a vector
#' @export
#'
#' @examples
#' norm1 (rnorm(10))


norm1 <- function(vect) {
  sum(abs(vect))
}

#roxygen2::roxygenise()
