#' GSVD
#'
#' @param Y A data frame or matrix with factors
#' @param R The number of components
#'
#' @return Returns the R first generalized singular triplets (U, V, d) of a data frame with factors
#' @export
#'
#' @examples
#' library(FactoMineR)
#' data(tea)
#' tea <- tea [, 1:8]
#' gsvd(tea, 3)

gsvd <- function (Y, X, R = 5, row.w = NULL) {
  #tableau disjonctif des donnees
  X <- as.matrix(X)

  # n <- ncol(X)- ncol(Y) #nb de vp max #ajouter une condition

  ##creation vecteurs poids et masses
  # if (is.null(row.w)) {
  #   r <- X %*% matrix(rep(1, dim(X)[2])) / sum(X)
  #   c <- t(X) %*% matrix(rep(1, dim(X)[1])) / sum(X)
  #   #centrer la matrice
  #   X <- X / sum(X) - r %*% t(c)
  # }else {
  #   mat <- as.matrix(X)*row.w/sum(X*row.w)
  #   # r <- row.w /sum(row.w)
  #   r <- as.matrix(rowSums(mat))
  #   c <- as.matrix(colSums(mat))
  #   #centrer la matrice
  #   X <- mat - r %*% t(c)
  #   }


  r <- X %*% matrix(rep(1, dim(X)[2])) / sum(X)
  c <- t(X) %*% matrix(rep(1, dim(X)[1])) / sum(X)

  # M <- diag(as.numeric(r)) #masses
  # W <- diag(as.numeric(c)) #weights

  # #centrer la matrice
  X <- X / sum(X) - r %*% t(c)

  #contrained matrix X'
  Xtilde <- diag((as.numeric(r) ^ (- 1 / 2))) %*% X %*% diag((as.numeric(c) ^ (- 1 / 2))) #(M %^% (-1 / 2)) %*% X %*% (W %^% (-1 / 2)) #constrained matrix

  res <- svd(Xtilde)

  Qtilde <- res$v
  D <- res$d^2
  Ptilde <- res$u

  Q <- diag(as.numeric(c) ^ (1 / 2)) %*% Qtilde#t(t(Qtilde) %*% (W %^% (-1 / 2)))
  P <- diag(as.numeric(r) ^ (1 / 2)) %*% Ptilde #t(t(Ptilde) %*% (M %^% (-1 / 2)))

  return (list(D = D, P = P, Q = Q, r = r, c = c, Xcent = X))
}
