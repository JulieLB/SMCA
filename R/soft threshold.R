#' Soft Threshold norme 1
#'
#' @param x a vector
#' @param alpha a numeric threshold
#'
#' @return Returns the proximal vector of x for a threshold alpha.
#' @export
#'
#' @examples
#' #the proximal vector of a random vector
#' soft_threshold(rnorm(10), 0.2)

soft_threshold <- function(x, alpha) {
  return(sign(x) * pmax(0, abs(x) - alpha))
  # ifelse(abs(x) < alpha, 0 , x - alpha * sign(x))
}
