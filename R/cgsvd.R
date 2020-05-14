#' CGSVD
#'
#' @param Y A data frame or matrix with factors
#' @param c1 The radius of l1 ball (constraint)
#' @param c2 The radius of l2 ball (constraint)
#' @param R The number of singular triplets we want to obtain
#' @param eps.pi Precision
#' @param itermax.pi The maximum number of iteration of power iteration
#' @param eps.pocs Precision
#' @param itermax.pocs The maximum number of iteration of projections
#' @param v.partition If true, creates a group constraint based on the ariables partition
#' @param init Method to initialize the singular vectors, by default "gsvd" (can also be "rand")
#' @param Gcol A group partition of the data columns (categories)
#' @param Grow A group partition of the data rows
#'
#' @return Returns the constrained singular triplets of factor data and the number of iteration used to obtain them.
#' @export
#'
#' @examples
#' library(FactoMineR)
#' data(tea)
#' tea <- tea [, 1:8]
#'
#' test <- cgsvd (tea, 1.5, 1.5, R = 5, init = "rand")

cgsvd <- function (Y,
                   X,
                   c1,
                   c2,
                   R = 2,
                   init = "svd",
                   eps.pi = 1e-12,
                   itermax.pi = 500,
                   eps.pocs = 1e-12,
                   itermax.pocs = 500,
                   v.partition = F,
                   Gcol = NULL,
                   Grow = NULL) {


  #tableau disjonctif des donnees
  X <- as.matrix(X)

  n <- ncol(X) - ncol(Y) #nb de vp que l'on cherche

  #partition
  if (v.partition) {Gcol <- partition_variables(Y)
  }else{ Gcol <- Gcol}

  ##creation vecteurs poids et masses
  r <- X %*% matrix(rep(1, dim(X)[2])) / sum(X)
  c <- (t(X) %*% matrix(rep(1, dim(X)[1])) / sum(X))

  M <- diag(as.numeric(r)) #masses
  W <- diag(as.numeric(c)) #weights

  #centrer la matrice
  X <- X / sum(X) - r %*% t(c)

  #contrained matrix X'
  Xtilde <- diag((diag(M) ^ (-1 / 2))) %*% X %*% diag((diag(W) ^ (-1 / 2))) #(M %^% (-1 / 2)) %*% X %*% (W %^% (-1 / 2)) #constrained matrix

  # X_csvd <- t(Xtilde) %*% Xtilde
  X_csvd <- Xtilde

  res <- csvd(X_csvd, R, c1 = c1, c2 = c2, init = init, eps.pi = eps.pi, itermax.pi = itermax.pi, eps.pocs = eps.pocs, itermax.pocs = itermax.pocs, Gcol = Gcol, Grow = Grow)

  Qtilde <- res$Q
  D <- res$D
  Ptilde <- res$P #data.matrix(Xtilde %*% Qtilde %*% diag(sqrt(D) ^ -1))

  Q <- t(t(Qtilde) %*% diag((diag(W) ^ (-1 / 2))))#t(t(Qtilde) %*% (W %^% (-1 / 2)))
  P <- t(t(Ptilde) %*% diag((diag(M) ^ (-1 / 2)))) #t(t(Ptilde) %*% (M %^% (-1 / 2)))

  return(list(D = D, P = P, Q = Q, iter = res$it, r = r, c = c))
}
