% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_SMCA.R
\name{plot_SMCA}
\alias{plot_SMCA}
\title{Sparse MCA plot}
\usage{
plot_SMCA(res, axes = c(1, 2), choix = "ind", aff.noms = F, habillage = NULL)
}
\arguments{
\item{res}{the result of SMCA or SPCA}

\item{axes}{The axes we want to show}

\item{choix}{With this options, we decide if we want the individual plot (choix = "ind") or the Categories plot (choix = "var")}

\item{aff.noms}{If TRUE (default), the points names are shown}

\item{habillage}{The number of the variable according to which we want to color the points. If NULL (default), all the points will be of the same color.}
}
\value{
Returns a map (ggplot2) of the individuals or the categories with options.
}
\description{
Sparse MCA plot
}
\examples{
data(cheese)
res <- SMCA(cheese, c1=10, c2 = 4, R = 10, meth = 'cgsvd')
plot.SMCA(res, choix = "ind", aff.noms = F, habillage = 3)
plot.SMCA(res, choix = "var", aff.noms = T)
}
