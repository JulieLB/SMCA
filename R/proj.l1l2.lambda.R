
#' Lambda argument of the projection into the intersection of l1 and l2 Balls
#'
#' @param x A vector
#' @param c A constraint (numeric)
#'
#' @return the lamdba argument of the proximal operator
#' @export
#'
#' @examples
#' proj.l1l2.lambda(rnorm(10), 1.5)

proj.l1l2.lambda <- function(x, c) {
  n <- length(x)
  if (c >= sqrt(n)) return(0)
  absx <- abs(x)
  nmax <- sum(absx == max(absx))
  # print(nmax)
  # print(sqrt(nmax))
  if (c <= sqrt(nmax)) {
    res <- rep(0, length(x))
    imax <- which(absx == max(absx))
    res[imax] <- x[imax]
    return(0)
  }

  #step 1 : valeur absolue dans l'ordre decroissant
  xtilde <- sort(abs(x), decreasing = T)

  # step2 : trouver i
  # on calcule psi pour chaque element de xtilde
  # m <- sapply(1:n, function(k) {
  #   psi(x, xtilde[k])
  # })
  m <- psi(x, xtilde)

  # on prend la plus grande valeur inferieure a c (contrainte)
  i <- max(which(m <= c))
  xtildei <- xtilde[i]
  # step 3 : calcul de delta
  term1 <- norm2(soft_threshold(xtilde, xtildei)) / i
  term2 <- (i - psi(x, xtildei) ^ 2) / (i - c ^ 2)
  d <- term1 * (c * sqrt(term2) - psi(x, xtildei))

  return(xtildei - d)
}
