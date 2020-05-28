#' Gives the Screeplot of a Sparce MCA or PCA
#'
#' @param res the result of SMCA or SPCA
#' @param aff.mean If TRUE, adds an horizontal line corresponding to the mean of the eigen values.
#' @param sort.comp If TRUE (default), the components are sorted by decreasing size.
#' @param title by default "Scree plot".
#' @param ncomp The number of components to show.
#'
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#'
#' @examples
#' data(cheese)
#' res.gpmd <- SMCA(cheese, n = 2, c1 = 3, c2 = 1.9, meth = "gpmd")
#' screeplot_new(res.gpmd)

screeplot_new <- function(res, aff.mean = F, sort.comp = T, title = "Scree plot", ncomp = "all") {

  if (ncomp == "all") { contrib <- data.frame(dim = res$eig$dim,
                                               eig = res$eig$percentageOfVariance)
  }else if (is.numeric(ncomp)) {contrib <- data.frame(dim = res$eig$dim[1:ncomp],
                                                      eig = res$eig$percentageOfVariance[1:ncomp])
  }else return("error")


  if (sort.comp) contrib$eig <- sort(contrib$eig, d = T)


  if(aff.mean) {
    aff.mean.plot <- geom_hline(yintercept = mean(res$eig$percentageOfVariance), color ="red")
  } else {
    aff.mean.plot <- NULL
  }

  screeplot <- ggplot2::ggplot(contrib, aes(dim, eig)) +
    geom_col(fill = 'lightblue') +
    geom_line() +
    geom_point() +
    aff.mean.plot +
    labs(title =  title,
                  x = "Components",
                  y = "Percentage of explained variance") +
    theme_light()

  return(screeplot)
}
