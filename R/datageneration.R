#' Generate Artificial Data by Simulation
#'
#' The `datageneration()` function generates artificial ring-network data
#' by simulation.
#' The function is used in the package vignette.
#'
#' @param n The sample size
#'
#' @returns A list containing the outcome vector, the treatment vector,
#' the instrumental vector, and the true instrumental exposure vector, and
#' the symmetric binary adjacency matrix.
#'
#' @examples
#' latenetwork::datageneration(n = 2000)
#'
#' @export
#'
datageneration <- function(n) {

  # Adjacency matrix of ring network
  ring <- igraph::make_ring(n,
                            directed = FALSE,
                            mutual = FALSE,
                            circular = TRUE)
  A <- igraph::as_adjacency_matrix(ring)
  A <- as.matrix(A)

  # Instrumental vector
  Z <- sample(x = 0:1,
              size = n,
              replace = TRUE,
              prob = c(0.5, 0.5))

  # Instrumental exposure
  IEM <- ifelse(A %*% Z > 0, 1, 0)

  # Treatment vector
  D <- ifelse(stats::rnorm(n, mean = -1.5, sd = 1) + Z + IEM >= 0, 1, 0)

  # Outcome vector
  Y <- stats::rnorm(n, mean = 1, sd = 1) + stats::runif(n, min = 1, max = 2) * D

  # Return
  return(list(Y = Y,
              D = D,
              Z = Z,
              IEM = IEM,
              A = A))

}
