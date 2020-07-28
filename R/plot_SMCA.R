#' Sparse MCA plot
#'
#' @param res the result of SMCA or SPCA
#' @param axes The axes we want to show
#' @param choix With this options, we decide if we want the individual plot (choix = "ind") or the Categories plot (choix = "var")
#' @param aff.noms If TRUE (default), the points names are shown
#' @param habillage The number of the variable according to which we want to color the individuals points. If NULL (default), all the points will be of the same color.
#' @param vect.habillage A vector containing the habillage you wish for you plot (it must be of the same length as the chosen plot number of points)
#' @param title.precision Character strings to add a precision to the title.
#' @param alpha Indicate the opacity of the factor map points (default 0.5).
#'
#' @return Returns a map (ggplot2) of the individuals or the categories with options.
#' @export
#'
#' @examples
#' data(cheese)
#' res <- SMCA(cheese, c1=10, c2 = 4, R = 10, meth = 'cgsvd')
#' plot.SMCA(res, choix = "ind", aff.noms = F, habillage = 3)
#' plot.SMCA(res, choix = "var", aff.noms = T)

plot_SMCA <- function (res, axes = c(1,2), choix = "ind", aff.noms = F, habillage = NULL, vect.habillage = NULL, title.precision = NULL, alpha = 0.5) {

  #est ce qu'il s'agit d'un plot sur les individus ou sur les variables
  if (choix == "ind") {
    coord <- res$ind$coord
    title <- "Individuals factor map"
    col <- "blue"
    shap <- 16
    size <- res$ind$contrib[, axes[1]] * res$eig$eigenvalue[axes[1]] + res$ind$contrib[, axes[2]] * res$eig$eigenvalue[axes[2]]
    ix_label <- which(res$ind$contrib[,axes[1]] > 0.3 | res$ind$contrib[,axes[2]] > 0.3)

  }else if (choix == "mod") {
    coord <- res$var$coord
    title <- "Categories factor map"
    col <- "red"
    shap <- 17
    size <- res$var$contrib[, axes[1]] * res$eig$eigenvalue[axes[1]] + res$var$contrib[, axes[2]] * res$eig$eigenvalue[axes[2]]
    ix_label <- which(res$var$contrib[, axes[1]] > mean(res$var$contrib[,axes[1]]) | res$var$contrib[, axes[2]] > mean(res$var$contrib[,axes[2]]))

    if (is.null(res$other$Gcol)) {Gcol <- partition_variables(res$other$Xinit)
    }else {Gcol <- res$other$Gcol}

    hab <- unlist(lapply(1:length(Gcol), function(i) {Gcol[[i]] <- rep(colnames(res$other$Xinit)[i], length(Gcol[[i]]))}))

  }else if (choix == "var") {
    coord <- res$var$eta2
    title <- "Variables factor map "
    size <- 3
    col <- "red"
    shap <- 18
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
    if (aff.noms) {aff.noms.plot <- ggrepel::geom_text_repel(size = 3, aes(color = factor(hab)), force = 2)}
  } else if (choix == "mod"){
    coord <- cbind (coord, hab = hab)
    aff.points <- geom_point(aes(size = size, colour = factor(hab)), shape = shap, alpha = alpha)
    if (aff.noms) {aff.noms.plot <- ggrepel::geom_text_repel(size = 3, aes(color = factor(hab)), force = 2)}
  } else if (!is.null(vect.habillage) && length(vect.habillage)==nrow(coord)){
    coord <- cbind (coord, hab = vect.habillage)
    aff.points <- geom_point(aes(size = size, colour = factor(hab)), shape = shap, alpha = alpha)
    if (aff.noms) {aff.noms.plot <- ggrepel::geom_text_repel(size = 3, aes(color = factor(hab)), force = 2)}
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
