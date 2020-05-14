
#' Principal Composant analysis
#'
#' @param X A data frame or matrix with numerics
#' @param n The number of components we want to obtain
#'
#' @return Returns the eigen values, factor scores ad csvd of a principal components analysis.
#' @export
#'
#' @examples
#' library(FactoMineR)
#' data(decathlon)
#' deca <- scale(decathlon[, -13], center=T)/sqrt(nrow(decathlon)-1)
#' PCA(deca, n = 4)

PCA <- function (X,  n = 5) {
  res <- svd (X)

  eig <- as.data.frame(cbind(dim = seq(from = 1, by = 1, length = ncol(X)),
                             eigenvalue = sort(res$d, d=T),
                             percentageOfVariance = sort(sapply (1:length(res$d), function(j){res$d[j]/sum(res$d)*100}), d = T),
                             cumulatedPercentageOfVariance = cumsum(sapply (1:length(res$d), function(j){res$d[j]/sum(res$d)*100}))
  )
  )
  F <- res$u %*% diag(sqrt(res$d))[,1:n]
  contrib <- (F^2) %*% diag(1/colSums(F^2)) * 100
  cos2 <- t(t(F^2)%*%diag(1/rowSums(F^2)))

  col <- paste("dim", seq(from = 1, by = 1, length = n))
  colnames(F) <- col; colnames(contrib) <- col; colnames (cos2) <- col
  rownames(F) <- rownames(contrib) <- rownames (cos2) <- rownames (X)
  ind <- list(coord = F, contrib = contrib, cos2 = cos2)


  G <- res$v %*% diag(sqrt(res$d))[,1:n]
  contrib <- (G^2) %*% diag(1/colSums(G^2)) * 100
  cos2 <- t(t(G^2)%*%diag(1/rowSums(G^2)))

  colnames(G) <- colnames(contrib) <- colnames (cos2) <- paste("dim", seq(from = 1, by = 1, length = n))
  rownames(G) <- rownames(contrib) <- rownames (cos2) <- colnames (X)
  var <- list(coord = G, contrib = contrib, cos2 = cos2)

  return(list (svd = res, eig = eig, var = var, ind = ind))
}
