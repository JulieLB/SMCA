#' Power iteration
#'
#' @param X A matrix of numerics
#' @param c1 the radius of l1 ball
#' @param c2 the radius of l2 ball
#' @param P0 Initial U values (random numeric vector or SVD first component)
#' @param Q0 Initial V values (random numeric vector or SVD first component)
#' @param P A matrix containing the previous singular values U
#' @param Q A matrix containing the previous singular values V
#' @param itermax The maximum number of iteration
#' @param eps Precision
#'
#' @return Returns the singular triplet of a matrix, given the already found singular triplets.
#' @export
#'
#' @examples
#' data(cheese)
#' X <- tab_disjonctif(cheese)
#' c1 <- 4
#' c2 <- 1.5
#' P0 <- normalize(rnorm(dim(X)[1]))
#' Q0 <- normalize(rnorm(dim(X)[2]))
#' PI (X, c1, c2, P0, Q0)

PI <- function(X,
               c1,
               c2,
               P0,
               Q0,
               P=NULL,
               Q=NULL,
               eps.pi = 1e-16,
               itermax.pi = 1000,
               eps.pocs = 1e-16,
               itermax.pocs = 1000,
               Gcol = NULL,
               Grow = NULL) {

  pnew <- pold <- as.vector(P0) #normalize(rnorm(nrow(X)))
  qnew <- qold <- as.vector(Q0) #normalize(rnorm(ncol(X)))

  for (j in 1:itermax.pi) {
    qnew <- proj.l1l2.orth(t(X) %*% pold, M = Q, c = c2, itermax.pocs = itermax.pocs, eps.pocs = eps.pocs, G = Gcol)$x
    pnew <- proj.l1l2.orth(X %*% qnew, M = P, c = c1, itermax.pocs = itermax.pocs, eps.pocs = eps.pocs, G = Grow)$x

    if ((norm2(pold - pnew) < eps.pi) && (norm2(qold - qnew) < eps.pi)) break

    pold <- pnew
    qold <- qnew
  }
  return( list(p=pnew, q=qnew, iter=j))
}

