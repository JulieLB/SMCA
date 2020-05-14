#' Soft Threshold norme de groupe
#'
#' @param x a vector
#' @param alpha a numeric threshold
#' @param G A group partition of the vector x
#'
#' @return Returns the group norm proximal vector of x for a threshold alpha.
#' @export
#'
#' @examples
#' #the proximal vector of a random vector
#' G <- list(1:4, 5, 6:8, 9:10)
#' x <- rnorm(10)
#' c <- 1.5
#' lambda <- proj.l1l2.lambda(x, c)
#' soft_threshold_normg(x, lambda, G)

soft_threshold_normg <- function(x, alpha, G) {

  #if the partition does not match to the vector, return the norm1 soft threshold
  if (length(x) != sum(sapply(1:length(G), function (i) {length(G[[i]])}))) {
      warning("Lengths of vector and group do not match!")
      return(soft_threshold(x, alpha))
    }

  prox <- rep(0, length(x))

  for (i in 1:length(G)) {
    iota <- G[[i]]
    xg <- normalize(x[iota]) * pmax(0, norm2(x[iota]) - alpha)
    prox[iota] <- xg
  }

  return(prox)
}
