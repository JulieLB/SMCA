#' Gives the Screeplots of multiple Sparce MCA or PCA
#'
#' @param res A data frame or matrix with the eigen values of different analysis as columns
#' @param title by default "Scree plot".
#' @param ncomp The number of components to show.
#'
#' @return A ggplot2 object
#' @export
#' @import ggplot2 and tidyverse
#'
#' @examples
screeplot_cum <- function(res, title = "Scree plot", ncomp = "all", names.res = NULL) {

  if(!is.data.frame(res)) return("The object res must be a data frame or matrix containing the eigenvalues.")

  if (ncomp == "all") { contrib <- cbind(dim = 1:nrow(res), res)
  }else if (is.numeric(ncomp)) {contrib <- cbind(dim = 1:ncomp, res[1:ncomp,])
  }else return("error")

  if (is.null(names.res)) {names <- sapply(1:ncol(res), function(i) {paste("res", i, sep = "")})
  }else{names <- names.res}

  colnames(contrib) <- c("dim", names)

  df <- contrib %>%
    select(dim, names) %>%
    gather(key = "variable", value = "value", -dim)

  screeplot <- ggplot2::ggplot(df, aes(x = dim, y = value)) +
    geom_line(aes(color = variable), size = 1) +
    geom_point(aes(color = variable)) +
    labs(title =  title,
         x = "Components",
         y = "Percentage of explained variance") +
    #scale_color_manual(values = c("darkred", "steelblue"))
    theme_light()

  return(screeplot)
}
