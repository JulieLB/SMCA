#' Orthonormal projection
#'
#' @param x A vector to project.
#' @param M A matrix. The vector is projected into the orthonormal space of M. #the matrix into which orthonormal space the vector is projected.
#'
#' @return The vector x projected into the orthonormal space of M.
#' @export
#'
#' @examples
#' x <- rnorm(10)
#' A <- matrix(data = c(1, 0.2, -0.3, 0.3, 0.5), nrow=10, ncol=3)
#' proj.orth(x, A)

# proj.orth <- function(x, M) {
proj.orth <- function(x, MtM) {
  x <- as.vector(x)
  # x <- x - M%*%t(M)%*%x
  x <- x - MtM %*% x
  return(as.vector(x))
}

