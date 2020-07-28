
#' Projection into intersection orthonormal matrix and l1 and l2 balls
#'
#' @param x A vector.
#' @param c A constraint (numeric).
#' @param itermax.pocs
#' @param eps.pocs
#' @param G A group partition of the vector x
#' @param M A matrix.
#'
#' @return The projected vector x into the intersection of the orthonormal space of M and the l1 and l2 Balls intersection.
#' @export
#'
#' @examples
#'

proj.l1l2.orth <- function(x, c, M = NULL, itermax.pocs = 1000, eps.pocs = 1e-16, G = NULL, orth.first = T) {
  if (is.null(M)) {
    return(list(x = proj.l1l2(x, c = c), it = 1))
  }
  xold <- xnew <- x
  MtM <- M %*% t(M) # ?crossprod

  if(orth.first){
    if (is.null(G)) {
      for (k in 1:itermax.pocs) {
        # xnew <- proj.l1l2(proj.orth(x = xold, M), c = c)
        xnew <- proj.orth(x = proj.l1l2(x = xold, c = c), MtM)
        if (norm2(xnew - xold) < eps.pocs) break
        xold <- xnew
      }
    }else{
      for (k in 1:itermax.pocs) {
        # xnew <- proj.lgl2(proj.orth(x = xold, M), c = c, G = G)
        xnew <- proj.orth(x = proj.lgl2(x = xold, c = c, G = G), MtM)
        if (norm2(xnew - xold) < eps.pocs) break
        xold <- xnew
      }
    }

  }else{
    if (is.null(G)) {
      for (k in 1:itermax.pocs) {
        # xnew <- proj.l1l2(proj.orth(x = xold, M), c = c)
        xnew <- proj.l1l2(proj.orth(x = xold, MtM), c = c)
        if (norm2(xnew - xold) < eps.pocs) break
        xold <- xnew
      }
    }else{
      for (k in 1:itermax.pocs) {
        # xnew <- proj.lgl2(proj.orth(x = xold, M), c = c, G = G)
        xnew <- proj.lgl2(proj.orth(x = xold, MtM), c = c, G = G)
        if (norm2(xnew - xold) < eps.pocs) break
        xold <- xnew
      }
    }
  }

  # if (is.null(G)) {
  #   for (k in 1:itermax.pocs) {
  #     # xnew <- proj.l1l2(proj.orth(x = xold, M), c = c)
  #     xnew <- proj.l1l2(proj.orth(x = xold, MtM), c = c)
  #     if (norm2(xnew - xold) < eps.pocs) break
  #     xold <- xnew
  #   }
  # }else{
  #   for (k in 1:itermax.pocs) {
  #     # xnew <- proj.lgl2(proj.orth(x = xold, M), c = c, G = G)
  #     xnew <- proj.lgl2(proj.orth(x = xold, MtM), c = c, G = G)
  #     if (norm2(xnew - xold) < eps.pocs) break
  #     xold <- xnew
  #   }
  # }



  return(list(x = xnew, it = k))
}
