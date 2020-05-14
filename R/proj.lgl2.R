#' The projection into the intersection of l12 and l2 Balls
#'
#' @param x A vector
#' @param c A constraint (numeric)
#' @param G A group partition of the vector x
#'
#' @return Returns the vector x projected into B12(c)interB2(1).
#' @export
#'
#' @examples
#' G <- list(1:3, 4:5)
#' x <- rnorm(5)
#' c <- 2*sqrt(5)/3
#' proj.l12l2(x, c, G)

proj.lgl2 <- function(x, c, G) {
  v <- v_subvector(x, G)
  # If c <= sqrt(nmax) then return...
  absv <- abs(v)
  nmax <- sum(absv == max(absv))

  if (c <= sqrt(nmax)) {
    gpmax <- which(absv == max(absv))
    x_prox <- rep (0, length(x))

    for (i in 1:length(gpmax)) {
      iota <- G[[i]]
      x_prox[iota] <- normalize(x[iota])
    }
  }else{
    lambda <- proj.l1l2.lambda(v, c)
    x_prox <- soft_threshold_normg(x, lambda, G)
  }

  return(normalize(x_prox))
}
