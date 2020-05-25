#' Sparse MCA plot
#'
#' @param res the result of SMCA or SPCA
#' @param axes The axes we want to show
#' @param choix With this options, we decide if we want the individual plot (choix = "ind") or the Categories plot (choix = "var")
#' @param aff.noms If TRUE (default), the points names are shown
#' @param habillage The number of the variable according to which we want to color the individuals points. If NULL (default), all the points will be of the same color.
#'
#' @return Returns a map (ggplot2) of the individuals or the categories with options.
#' @export
#'
#' @examples
#' data(cheese)
#' res <- SMCA(cheese, c1=10, c2 = 4, R = 10, meth = 'cgsvd')
#' plot.SMCA(res, choix = "ind", aff.noms = F, habillage = 3)
#' plot.SMCA(res, choix = "var", aff.noms = T)

plot_SMCA <- function (res, axes = c(1,2), choix = "ind", aff.noms = F, habillage = NULL, title.precision = NULL) {

  #est ce qu'il s'agit d'un plot sur les individus ou sur les variables
  if (choix == "ind") {
    coord <- res$ind$coord
    title <- "Individuals factor map"
    col <- "blue"
    shap <- 16
    alpha <- 0.5
    size <- res$ind$contrib[, axes[1]] * res$eig$eigenvalue[axes[1]] + res$ind$contrib[, axes[2]] * res$eig$eigenvalue[axes[2]]
    ix_label <- which(res$ind$contrib[,axes[1]] > 0.3 | res$ind$contrib[,axes[2]] > 0.3)

  }else if (choix == "mod") {
    coord <- res$var$coord
    title <- "Categories factor map"
    col <- "red"
    shap <- 17
    alpha <- 0.5
    size <- res$var$contrib[, axes[1]] * res$eig$eigenvalue[axes[1]] + res$var$contrib[, axes[2]] * res$eig$eigenvalue[axes[2]]
    ix_label <- which(res$var$contrib[,axes[1]] > 0.3 | res$var$contrib[,axes[2]] > 0.3)

  }else if (choix == "var") {
    coord <- res$var$eta2
    title <- "Variables factor map "
    size <- 3
    col <- "red"
    shap <- 18
    alpha <- 0.5
    ix_label <- which(res$var$eta2[,axes[1]] > 0.001 | res$var$eta2[,axes[2]])
  }

  #add a precision to the title ?
  if (is.character(title.precision)) {
    title <- paste (title, title.precision, sep = " ")
  }


  contrib <- res$eig$percentageOfVariance
  comp.names <- sapply(1:dim(coord)[2], function (i) {paste("Comp", i, "  (", round(contrib[i], 2), "%", ")", sep="")})

  coord <- as.data.frame(cbind(x = coord[, axes[1]], y = coord[, axes[2]], names = NA))

  #est ce qu'on affiche le nom des points
  aff.noms.plot <- NULL
  if (aff.noms==T) {
    coord$names[ix_label] <- rownames(coord) [ix_label]
    aff.noms.plot <- ggrepel::geom_text_repel(size = 3, color = col, force = 2) #, position=position_jitter(width=0.2,height=0.2)) # vjust = "outward", hjust = "outward")
  }

  #habillage
  aff.points <- geom_point(aes(size = size), shape = shap, colour = col, alpha = alpha)
  if (!is.null(habillage) && (habillage %in% 1:ncol(res$other$Xinit))) {
    coord <- cbind (coord, hab = res$other$Xinit[,habillage])
    aff.points <- geom_point(aes(size = size, colour = factor(hab)), shape = shap, alpha = alpha)
  }

  coord$x <- as.numeric(as.vector(coord$x))
  coord$y <- as.numeric(as.vector(coord$y))


  ggplot(coord, aes(x, y, label = names)) +
    aff.points +
    aff.noms.plot +
    labs(title = title, x = comp.names[axes[1]], y = comp.names[axes[2]]) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # expand_limits(x = c(min(coord$x)-1, max(coord$x)+1), y = c(min(coord$y)-1, max(coord$y)+1))+
    theme_light()

}
