
#' Multivariate Correspondance Analysis
#'
#' @param Y A data frame or matrix with factors
#' @param n The number of components
#'
#' @return Returns the R first factor scores of Y
#' @export
#'
#' @examples
#' data(cheese)
#' MCA(cheese, 3)

MCA <- function(Y, n = 5) {

  if(sum(is.na(Y))>0) stop("Error. Missing values.")

  X <- tab_disjonctif(Y)

  if (n > min(ncol(X), nrow(X))) n <- min(ncol(X), nrow(X))

  res <- gsvd(Y = Y, X = X, R = n)

  Gcol <- partition_variables(Y)
  Itot <- 1/length(Gcol)*sum(sapply(1:length(Gcol), function(i) {length(Gcol[[i]])-1}))

  eig <- as.data.frame(cbind(dim = seq(from = 1, by = 1, length = n),
                             eigenvalue = res$D[1:n],
                             percentageOfVariance = res$D[1:n] / Itot * 100,
                             cumulatedPercentageOfVariance = cumsum(res$D[1:n] / Itot * 100)
  )
  )

  #individus
  F <- diag(as.numeric(res$r)^(-1)) %*% res$P %*% diag(sqrt(res$D))[,1:n] #diag(as.numeric(res$r)^(-1/2)) %*% res$P %*% diag(sqrt(res$D))[,1:n]
  F2 <- F^2
  contrib <- 100/sum(X) * diag(rowSums(X)) %*% F2 %*% diag(1/res$D[1:n]) #(F2) %*% diag(1/colSums(F2)) * 100 #
  cos2 <- t(t(F2)%*%diag(1/rowSums(F2)))

  col <- paste("dim", seq(from = 1, by = 1, length = n))
  colnames(F) <- colnames(contrib) <- colnames (cos2) <- col
  rownames(F) <- rownames(contrib) <- rownames (cos2) <- rownames (Y)
  ind <- list(coord = F, contrib = contrib, cos2 = cos2)

  #categories
  G <- diag(as.numeric(res$c)^(-1)) %*% res$Q %*% diag(sqrt(res$D))[,1:n] #diag(as.numeric(res$c)^(-1/2)) %*% res$Q %*% diag(sqrt(res$D))[,1:n]
  G2 <- G^2
  contrib <- 1/sum(X) * diag(colSums(X)) %*% G2 #contribution absolue
  cos2 <- t(t(G2)%*%diag(1/rowSums(G2)))

  partition <- partition_variables(Y)

  eta2 <- c()
  for (j in 1:n) {
    eta2 <- cbind(eta2, sapply(1:ncol(Y), function(i) {sum(contrib[partition[[i]],j])}))
  }

  colnames(G) <- colnames(contrib) <- colnames (cos2) <- colnames(eta2) <- col
  rownames(G) <- rownames(contrib) <- rownames (cos2) <- colnames (X)
  rownames(eta2) <- colnames(Y)
  var <- list(coord = G, contrib = 100 * contrib %*% diag(1/res$D[1:n]), cos2 = cos2, eta2 = 10 * eta2)


  return(list (gsvd = res,
               eig = eig,
               var = var,
               ind = ind,
               other =list(Xinit = Y,
                           Xdisj = X)))
}
