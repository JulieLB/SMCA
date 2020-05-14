
#' Power iteration first component
#'
#' @param X A matrix of numerics
#' @param c1 the radius of l1 ball
#' @param c2 the radius of l2 ball
#' @param U0 Initial U values (random numeric vector or SVD first component)
#' @param V0 Initial V values (random numeric vector or SVD first component)
#' @param itermax The maximum number of iteration
#' @param eps Precision
#'
#' @return
#' @export
#'
#' @examples
#'

PI.1 <- function(X, c1, c2, U0, V0, itermax=1000, eps=1e-16) {
  unew <- uold <- as.vector(U0) #normalize(rnorm(nrow(X)))
  vnew <- vold <- as.vector(V0) #normalize(rnorm(ncol(X)))

  for (j in 1:itermax) {

    unew <- normalize( proj.l1l2(X %*% vold, c1) )
    vnew <- normalize( proj.l1l2(t(X) %*% unew, c2) )

    if ((norm2(uold - unew) < eps) && (norm2(vold - vnew) < eps)) break
    uold <- unew
    vold <- vnew
  }
  return( list(u=unew, v=vnew, iter=j))
}
