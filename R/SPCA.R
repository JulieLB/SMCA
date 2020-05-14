#' Sparse Principal Components Analysis
#'
#' @param X A data frame or matrix with numerics
#' @param c1 The radius of l1 ball (constraint). It can take values from 1 to the square root of the number of row. The closer to 1 it is, the more sparsified are the observations.
#' @param c2 The radius of l2 ball (constraint). It can take values from 1 to the square root of the number of columns. The closer to 1 it is, the more sparsified are the variables.
#' @param n The number of components we want to obtain
#'
#' @return Returns the eigen values, factor scores ad csvd of a sparse principal components analysis.
#' @export
#'
#' @examples
#' library(FactoMineR)
#' data(decathlon)
#' deca <- scale(decathlon[, -13], center=T)/sqrt(nrow(decathlon)-1)
#' SPCA(deca, c1 = sqrt(nrow(deca))/2, c2 = sqrt(ncol(deca))/2, n = 4)

SPCA <- function (X, c1, c2, n) {
  res <- csvd (X = X, c1 = c1, c2 = c2, R = n)

  eig <- as.data.frame(cbind(dim = seq(from = 1, by = 1, length = n),
                             eigenvalue = sort(res$D, d=T),
                             percentageOfVariance = sort(sapply (1:length(res$D), function(j){res$D[j]/sum(res$D)*100}), d = T),
                             cumulatedPercentageOfVariance = cumsum(sapply (1:length(res$D), function(j){res$D[j]/sum(res$D)*100}))
  )
  )

  F <- res$P %*% diag(sqrt(res$D))[,1:n]
  contrib <- (F^2) %*% diag(1/colSums(F^2)) * 100
  cos2 <- t(t(F^2)%*%diag(1/rowSums(F^2)))

  col <- paste("dim", seq(from = 1, by = 1, length = n))
  colnames(F) <- colnames(contrib) <- colnames (cos2) <- col
  rownames(F) <- rownames(contrib) <- rownames (cos2) <- rownames (X)
  ind <- list(coord = F, contrib = contrib, cos2 = cos2)


  G <- res$Q %*% diag(sqrt(res$D))[,1:n]
  contrib <- (G^2) %*% diag(1/colSums(G^2)) * 100
  cos2 <- t(t(G^2)%*%diag(1/rowSums(G^2)))

  colnames(G) <- colnames(contrib) <- colnames (cos2) <- col
  rownames(G) <- rownames(contrib) <- rownames (cos2) <- colnames (X)
  var <- list(coord = G, contrib = contrib, cos2 = cos2)

  return(list (csvd = res, eig = eig, var = var, ind = ind))
}
