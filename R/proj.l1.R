#' Projection into l1 Ball
#'
#' @param x A vector.
#' @param c A constraint (numeric).
#'
#' @return Returns the projection of vector x into B1(c).
#' @export
#'
#' @examples
#' proj.l1(rnorm(10), 1.5)

proj.l1 <- function (x, c) {
  xtilde <- sort(abs(x), decreasing=T)
  n <- length(x)
  m <- sapply(1:n, function (k) {norm1(soft_threshold(x, xtilde[k]))})
  i <- max(which(m<=c))
  d <- (c-norm1(soft_threshold(x, xtilde[i])))/i
  return (soft_threshold(x, xtilde[i]-d))
}
