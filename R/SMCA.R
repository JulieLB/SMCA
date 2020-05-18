#' Sparse Multivariate Correspondance Analysis
#'
#' @param Y A data frame or matrix with factors
#' @param c1 The radius of l1 ball (constraint). It can take values from 1 to the square root of the number of row. The closer to 1 it is, the more sparsified are the observations.
#' @param c2 The radius of l2 ball (constraint). It can take values from 1 to the square root of the number of columns. The closer to 1 it is, the more sparsified are the variables.
#' @param n The number of components we want to obtain
#' @param meth Method to create the sparse singular vectors, by default "cgsvd" (can also be "gpmd")
#' @param init Method to initialize the singular vectors, by default "rand" (can also be "svd")
#' @param v.partition If true, creates a group constraint based on the ariables partition (by default False)
#' @param Gcol A group partition of the data columns (categories)
#' @param Grow A group partition of the data row
#'
#' @return Returns the sparse MCA for the given constraints.
#' @export
#'
#' @examples
#'data(cheese)
#'SMCA(cheese, c1 = sqrt(nrow(cheese))/2, c2 = sqrt(ncol(cheese))/2, n = 4)

SMCA <- function(Y, c1, c2, n = 5, meth ='cgsvd', init = "rand", v.partition = F, Grow = NULL, Gcol = NULL) {
  X <- tab_disjonctif(Y)

  if (meth =='cgsvd') {
    res <- cgsvd(Y = Y, X = X, c1 = c1, c2 = c2, R = n, init = init, v.partition = v.partition, Grow = Grow, Gcol = Gcol)
  }else if (meth =='gpmd') {
    res <- gpmd(Y = Y, X = X, K = n, sumabsu = c1, sumabsv = c2)
  }

  eig <- as.data.frame(cbind(dim = seq(from = 1, by = 1, length = length(res$D)),
                             eigenvalue = res$D,
                             percentageOfVariance = sapply (1:length(res$D), function(j){res$D[j]/sum(res$D)*100}),
                             cumulatedPercentageOfVariance = cumsum(sapply (1:length(res$D), function(j){res$D[j]/sum(res$D)*100}))
                             )
                       )
  F <- res$P %*% diag(sqrt(res$D))[,1:n] #diag(as.numeric(res$r)^(-1/2)) %*% res$P %*% diag(sqrt(res$D))[,1:n]
  contrib <- (F^2) %*% diag(1/colSums(F^2))*100
  cos2 <- t(t(F^2)%*%diag(1/rowSums(F^2)))

  col <- paste("dim", seq(from = 1, by = 1, length = n))
  colnames(F) <- colnames(contrib) <- colnames (cos2) <- col
  rownames(F) <- rownames(contrib) <- rownames (cos2) <- rownames (Y)
  ind <- list(coord = F, contrib = contrib, cos2 = cos2)


  G <- res$Q %*% diag(sqrt(res$D))[,1:n] #diag(as.numeric(res$c)^(-1/2)) %*% res$Q %*% diag(sqrt(res$D))[,1:n]
  contrib <- (G^2) %*% diag(1/colSums(G^2))*100
  cos2 <- t(t(G^2)%*%diag(1/rowSums(G^2)))

  partition <- partition_variables(Y)
  eta2 <- c()
  for (j in 1:n) {
    eta2 <- cbind(eta2, sapply(1:ncol(Y), function(i) {sum((G[partition[[i]],j])^2)/sum(G[,j]^2)}))
  }

  colnames(G) <- colnames(contrib) <- colnames (cos2) <- colnames(eta2) <- col
  rownames(G) <- rownames(contrib) <- rownames (cos2) <- colnames (X)
  rownames(eta2) <- colnames(Y)
  var <- list(coord = G, contrib = contrib, cos2 = cos2, eta2 = eta2)



  return(list (cgsvd = res,
               eig = eig,
               var = var,
               ind = ind,
               other =list(Xinit = Y,
                           Xdisj = X)))
}
