#' Generalized PMD
#'
#' @param Y A data frame or matrix with factors
#' @param sumabsu The constraint on the rows (numeric)
#' @param sumabsv The constraint on the columns (numeric)
#' @param K The number of singular triplets we want to obtain
#'
#' @return Returns the constrained singular triplets of factor data with the penalised matrix decomposition from the package PMA.
#' @export
#'
#' @examples
#' library(PMD)
#' library(FactoMineR)
#' data(tea)
#' tea <- tea [, 1:8]
#'
#' test <- gpmd (tea, 1.5, 1.5, K = 5)

gpmd <- function (Y, X, sumabsu, sumabsv, K = 3) {
  X <- as.matrix(X)

  K <- ncol(X) - ncol(Y) #nb de vp que l'on cherche

  ##creation vecteurs poids et masses
  r <- X %*% matrix(rep(1, dim(X)[2])) / sum(X)
  c <- (t(X) %*% matrix(rep(1, dim(X)[1])) / sum(X))

  M <- diag(as.numeric(r)) #masses
  W <- diag(as.numeric(c)) #weights

  #centrer la matrice
  X <- X / sum(X) - r %*% t(c)

  #contrained matrix X'
  Xtilde <- diag((diag(M) ^ (-1 / 2))) %*% X %*% diag((diag(W) ^ (-1 / 2))) #(M %^% (-1 / 2)) %*% X %*% (W %^% (-1 / 2)) #constrained matrix

  # X_pmd <- t(Xtilde) %*% Xtilde
  X_pmd <- Xtilde

  res <- PMD(X_pmd, K = K, sumabsu = sumabsu, sumabsv = sumabsv)

  Qtilde <- res$v
  D <- res$d
  Ptilde <- res$u # data.matrix(Xtilde %*% Qtilde %*% diag(sqrt(D) ^ -1))

  Q <- t(t(Qtilde) %*% diag((diag(W) ^ (-1 / 2))))#t(t(Qtilde) %*% (W %^% (-1 / 2)))
  P <- t(t(Ptilde) %*% diag((diag(M) ^ (-1 / 2)))) #t(t(Ptilde) %*% (M %^% (-1 / 2)))

  return(list(D = D, P = P, Q = Q, M = M, W = W))

}
