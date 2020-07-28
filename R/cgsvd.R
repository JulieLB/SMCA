#' Constrained Generalized Singular Value Decomposition
#'
#' @param Y A data frame or matrix containing factors
#' @param X ???
#' @param c1 radius of the l1 ball (left)
#' @param c2 radius of the l2 ball (right)
#' @param R number of singular triplets we want to obtain
#' @param eps.pi Precision
#' @param itermax.pi The maximum number of iteration of power iteration
#' @param eps.pocs Precision
#' @param itermax.pocs The maximum number of iteration of projections
#' @param v.partition If true, creates a group constraint based on the ariables partition
#' @param init Method to initialize the singular vectors, by default "gsvd" (can also be "rand")
#' @param Gcol A group partition of the data columns (categories)
#' @param Grow A group partition of the data rows
#' @param row.w
#' @param double.centering
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
                   Grow = NULL,
                   row.w = NULL,
                   double.centering = T) {


  #tableau disjonctif des donnees
  X <- as.matrix(X)

  n <- ncol(X) - ncol(Y) # nb de vp que l'on cherche

  # partition
  if (v.partition) {
    Gcol <- partition_variables(Y)
  } else {
    Gcol <- Gcol
  }

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
  # }

  r <- X %*% matrix(rep(1, dim(X)[2])) / sum(X)
  c <- (t(X) %*% matrix(rep(1, dim(X)[1])) / sum(X))

  # #centrer la matrice
  X <- X / sum(X) - r %*% t(c)

  # M <- diag(as.numeric(r)) #masses
  # W <- diag(as.numeric(c)) #weights



  #contrained matrix X'
  Xtilde <- diag(as.numeric(r) ^ (-1 / 2)) %*% X %*% diag(as.numeric(c) ^ (-1 / 2)) #(M %^% (-1 / 2)) %*% X %*% (W %^% (-1 / 2)) #constrained matrix

  # X_csvd <- t(Xtilde) %*% Xtilde
  X_csvd <- Xtilde

  res <- csvd(X_csvd, R, c1 = c1, c2 = c2, r = r, c = c, init = init, eps.pi = eps.pi, itermax.pi = itermax.pi, eps.pocs = eps.pocs, itermax.pocs = itermax.pocs, Gcol = Gcol, Grow = Grow, double.centering = double.centering)
  # res <- csvd(X_csvd, R, c1 = c1, c2 = c2, r = NULL, c = NULL, init = init, eps.pi = eps.pi, itermax.pi = itermax.pi, eps.pocs = eps.pocs, itermax.pocs = itermax.pocs, Gcol = Gcol, Grow = Grow)


  Qtilde <- res$Q
  D <- res$D^2
  Ptilde <- res$P #data.matrix(Xtilde %*% Qtilde %*% diag(sqrt(D) ^ -1))

  Q <- diag(as.numeric(c) ^ (1 / 2)) %*% Qtilde #t(t(Qtilde) %*% diag((diag(W) ^ (-1 / 2))) ) #t(t(Qtilde) %*% (W %^% (-1 / 2)))
  P <- diag(as.numeric(r) ^ (1 / 2)) %*% Ptilde #t(t(Ptilde) %*% diag((diag(M) ^ (-1 / 2))) ) #t(t(Ptilde) %*% (M %^% (-1 / 2)))

  return(list(D = D, P = P, Q = Q, iter = res$it, r = r, c = c))
}
