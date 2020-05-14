#' Intermediate Function of proj.l1l2
#'
#' @param x A vector
#' @param lambda Argument
#'
#' @return
#' @export
#'
#' @examples
#'
psi <- Vectorize(function(x, lambda) {
  # x <- sort(abs(x), decreasing = T)
  x_tilde <- soft_threshold(abs(x), lambda)
  norm1x <- norm1(x_tilde)
  if (norm1x < 1e-16) return(1)
  else return(norm1x / norm2(x_tilde))
}, vectorize.args = "lambda")
