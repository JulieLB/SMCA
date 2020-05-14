#' Normalize a vector
#'
#' @param x a vector
#'
#' @return Return the normalization of a vector (the vector divided by its l2 norm)
#' @export
#'
#' @examples normalize (rnorm(10))
#'

normalize <- function(x) {
  normx <- norm2(x)
  if (normx < 1e-16)
    return(0 * x)
  else
    x / norm2(x)
}
