#' CSVD
#'
#' @param X A matrix of numerics
#' @param c1 The radius of l1 ball (constraint)
#' @param c2 The radius of l2 ball (constraint)
#' @param init Method to initialize the singular vectors, by default "svd" (can also be "rand")
#' @param R The number of singular triplets we want to obtain
#' @param eps.pi Precision
#' @param itermax.pi The maximum number of iteration of power iteration
#' @param eps.pocs Precision
#' @param itermax.pocs The maximum number of iteration of projections
#' @param Gcol A group partition of the data rows
#' @param Grow A group partition of the data columns
#' @param r The row weight (default = NULL). Add r if you want orthogonality of the left singular vectors with the row weight.
#' @param c The column weight (default = NULL). Add c if you want orthogonality of the right singular vectors with the column weight.
#'
#' @return Returns the constrained singular triplets of numeric data and the number of iteration used to obtain them
#' @export
#'
#' @examples
#'library(FactoMineR)
#'data("decathlon")
#'deca <- as.matrix(decathlon[,-13])
#'c1 <- sqrt(nrow(deca))/2
#'c2 <- sqrt(ncol(deca))/2
#'csvd2(deca, c1, c2, R=2)


csvd <- function (X,
                  c1,
                  c2,
                  r = NULL,
                  c = NULL,
                  init = "svd",
                  R = 3,
                  eps.pi = 1e-12,
                  itermax.pi = 1000,
                  eps.pocs = 1e-12,
                  itermax.pocs = 1000,
                  Gcol = NULL,
                  Grow = NULL) {

  # R <- min(nrow(X), ncol(X))

  if (is.null(r) | is.null(c)) {
    P <- matrix(0, nrow = nrow(X), ncol = R+1)
    Q <- matrix(0, nrow = ncol(X), ncol = R+1)
  }else if (length(r)==nrow(X) & length(c)==ncol(X)){
    P <- cbind(r, matrix(0, nrow = nrow(X), ncol = R))
    Q <- cbind(c, matrix(0, nrow = ncol(X), ncol = R))
  }else stop("Error. Wrong attribution of r or c.")

  iter <- c()

  if (init=="svd") {
    res <- svd(X)
    P0 <- res$u
    Q0 <- res$v
  }else {
    P0 <- sapply(1:R, function (i) normalize(rnorm(nrow(X))))
    Q0 <- sapply(1:R, function (i) normalize(rnorm(ncol(X))))
  }

  res1 <- PI(X, c1[1], c2[1], P0 = P0[, 1], Q0 = Q0[, 1], P = P, Q = Q, eps.pi = eps.pi, itermax.pi = itermax.pi, eps.pocs = eps.pocs, itermax.pocs = itermax.pocs, Gcol = Gcol, Grow = Grow)
  P[, 2] <- res1$p
  Q[, 2] <- res1$q

  iter <- c(iter, res1$iter)

  if (R!=1) {
    for (i in 2:R) {
      res2 <- PI(X, c1[i], c2[i], P0[,i], Q0[,i], P = P[,1:i], Q = Q[,1:i], eps.pi = eps.pi, itermax.pi = itermax.pi, eps.pocs = eps.pocs, itermax.pocs = itermax.pocs, Gcol = Gcol, Grow = Grow)
      P[,i+1] <- res2$p
      Q[,i+1] <- res2$q
      iter <- c(iter, res2$iter)
    }
  }
  P <- P[, -1]
  Q <- Q[, -1]
  d <- drop(crossprod(P, X %*% Q))
  if (length(d) == 1) {
    D <- d
  } else {
    D <- diag(d)
  }
  return(list(P = P, Q = Q, D = D, iter = iter))
}

