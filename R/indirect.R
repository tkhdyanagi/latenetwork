#' Inference on Average Indirect Effect Parameters
#'
#' Inference on the average indirect effect of the IV on the outcome,
#' that on the treatment receipt, and the local average indirect effect
#' in the presence of network spillover of unknown form
#'
#' @details
#' The `indirect()` function estimates the average indirect effect of the IV
#' on the outcome, that on the treatment receipt, and
#' the local average indirect effect via inverse probability weighting
#' in the approximate neighborhood interference framework.
#' The function also computes the standard errors and the confidence intervals
#' for the target parameters based on the network HAC estimation and
#' the wild bootstrap.
#' For more details, see Hoshino and Yanagi (2023).
#' The lengths of `Y`, `D`, `Z`, `S` and
#' of the row and column of `A` must be the same.
#' `K` must be a positive integer.
#' `bw` must be `NULL` or a non-negative number.
#' `B` must be `NULL` or a positive number.
#' `alp` must be a positive number between 0 and 0.5.
#'
#' @param Y An n-dimensional outcome vector
#' @param D An n-dimensional binary treatment vector
#' @param Z An n-dimensional binary instrumental vector
#' @param S An n-dimensional logical vector to indicate whether each unit
#' belongs to the sub-population S
#' @param A An n times n symmetric binary adjacency matrix
#' @param K A scalar to indicate the range of neighborhood
#' used for constructing the interference set.
#' Default is 1.
#' @param bw A scalar of the bandwidth used for the HAC estimation and
#' the wild bootstrap.
#' If `bw = NULL`, the rule-of-thumb bandwidth proposed by Leung (2022) is used.
#' Default is NULL.
#' @param B The number of bootstrap repetitions.
#' If `B = NULL`, the wild bootstrap is skipped.
#' Default is NULL.
#' @param alp The significance level. Default is 0.05.
#'
#' @returns A data.frame containing the following elements:
#' \item{est}{The parameter estimate}
#' \item{HAC_SE}{The standard error computed by the network HAC estimation}
#' \item{HAC_CI_L}{The lower bound of the confidence interval computed by
#' the network HAC estimation}
#' \item{HAC_CI_U}{The upper bound of the confidence interval computed by
#' the network HAC estimation}
#' \item{wild_SE}{The standard error computed by the wild bootstrap}
#' \item{wild_CI_L}{The lower bound of the confidence interval computed by
#' the wild bootstrap}
#' \item{wild_CI_U}{The upper bound of the confidence interval computed by
#' the wild bootstrap}
#' \item{bw}{The bandwidth used for the HAC estimation
#' and the wild bootstrap}
#' \item{size}{The size of the subpopulation S}
#'
#' @examples
#' # Generate artificial data
#' set.seed(1)
#' n <- 2000
#' data <- latenetwork::datageneration(n = n)
#'
#' # Arguments
#' Y   <- data$Y
#' D   <- data$D
#' Z   <- data$Z
#' S   <- rep(TRUE, n)
#' A   <- data$A
#' K   <- 1
#' bw  <- NULL
#' B   <- NULL
#' alp <- 0.05
#'
#' # Estimation
#' latenetwork::indirect(Y = Y,
#'                       D = D,
#'                       Z = Z,
#'                       S = S,
#'                       A = A,
#'                       K = K,
#'                       bw = bw,
#'                       B = B,
#'                       alp = alp)
#'
#' @references Hoshino, T., & Yanagi, T. (2023).
#' Causal inference with noncompliance and unknown interference.
#' arXiv preprint arXiv:2108.07455.
#'
#' Leung, M.P. (2022).
#' Causal inference under approximate neighborhood interference.
#' Econometrica, 90(1), pp.267-293.
#'
#' @export
#'
indirect <- function(Y,
                     D,
                     Z,
                     S,
                     A,
                     K = 1,
                     bw = NULL,
                     B = NULL,
                     alp = 0.05) {

  # Error handling -------------------------------------------------------------

  error <- errorhandling(Y = Y,
                         D = D,
                         Z = Z,
                         IEM = NULL,
                         S = S,
                         A = A,
                         bw = bw,
                         B = B,
                         alp = alp)

  # Variable definitions -------------------------------------------------------

  # Distance matrix
  graph0 <- igraph::graph_from_adjacency_matrix(adjmatrix = A,
                                                mode = "undirected")

  distance0 <- igraph::distances(graph = graph0,
                                 v  = igraph::V(graph0),
                                 to = igraph::V(graph0),
                                 algorithm = "dijkstra")

  # Variable definitions
  Y0 <- Y
  D0 <- D
  Z0 <- Z

  Y <- Y[S]
  D <- D[S]
  Z <- Z[S]
  A <- A[S, S]
  graph <- igraph::graph_from_adjacency_matrix(adjmatrix = A,
                                               mode = "undirected")
  distance <- distance0[S, S]

  # Size of the sub-population S
  size <- sum(S)

  # List of neighbors
  neighbor_list <- list()
  for (i in which(S)) {
    neighbor <- which(distance0[i, ] >= 1 & distance0[i, ] <= K)
    neighbor_list <- append(neighbor_list, list(neighbor))
  }

  # Bandwidth -----------------------------------------------------------------

  if (is.null(bw)) {

    # Find largest connected component
    connect <- igraph::components(graph = graph)
    connect_member <-
      which(connect$membership == statip::mfv(connect$membership))

    # Average path length
    average_path_length <-
      sum(distance[connect_member, connect_member]) /
      (length(connect_member) * (length(connect_member) - 1))

    # Average degree
    average_degree <- sum(A) / size

    # Rule-of-thumb bandwidth proposed by Leung (2022)
    if (average_degree > 1) {

      bw0 <- ifelse(average_path_length < 2 * log(size) / log(average_degree),
                    0.5 * average_path_length,
                    average_path_length^(1/3))

      bw <- round(max(bw0, 2 * K))

    } else {

      bw <- 2 * K

    }

  }

  # Causal inference -----------------------------------------------------------

  # Function to compute the generalized propensity score
  # Input "z": A value of IV
  # Output: A value of GPS
  GPS <- function(z) {
    mean(Z == z)
  }

  # Function to compute the conditional outcome mean
  # Input "z": A value of IV
  # Output: A value of the conditional outcome mean
  mu_Y <- function(z) {
    summand <- rep(0, size)
    for (i in 1:size) {
      neighbor <- neighbor_list[[i]]
      if (length(neighbor) > 0) {
        summand[i] <- sum(Y0[neighbor]) * (Z[i] == z) / GPS(z = z)
      }
    }
    return(mean(summand))
  }

  # Function to compute the conditional treatment mean for AIED
  # Input "z": A value of IV
  # Output: A value of the conditional treatment mean for AIED
  mu_D_AIED <- function(z) {
    summand <- rep(0, size)
    for (i in 1:size) {
      neighbor <- neighbor_list[[i]]
      if (length(neighbor) > 0) {
        summand[i] <- sum(D0[neighbor]) * (Z[i] == z) / GPS(z = z)
      }
    }
    return(mean(summand))
  }

  # Function to compute the conditional treatment mean for ADED
  # Input "z": A value of IV
  # Output: A value of the conditional treatment mean for ADED
  mu_D_ADED <- function(z) {
    mean(D * (Z == z)) / GPS(z = z)
  }

  # Function to compute AIEY
  # Input: Nothing
  # Output: A value of AIEY
  AIEY <- function() {
    mu_Y(z = 1) - mu_Y(z = 0)
  }

  # Function to compute AIED
  # Input: Nothing
  # Output: A value of AIED
  AIED <- function() {
    mu_D_AIED(z = 1) - mu_D_AIED(z = 0)
  }

  # Function to compute ADED
  # Input: Nothing
  # Output: A value of ADED
  ADED <- function() {
    mu_D_ADED(z = 1) - mu_D_ADED(z = 0)
  }

  # Function to compute LAIE
  # Input: Nothing
  # Output: A value of LAIE
  LAIE <- function() {
    AIEY() / ADED()
  }

  # Estimates
  est_AIEY <- AIEY()
  est_AIED <- AIED()
  est_ADED <- ADED()
  est_LAIE <- LAIE()

  # Variable definitions
  W_Z1 <- (Z == 1)
  W_Z0 <- (Z == 0)

  W_Y <- W_D_AIED <- W_D_ADED <- rep(0, size)
  for (i in 1:size) {

    neighbor <- neighbor_list[[i]]

    multiplier <- (Z[i] == 1) / GPS(z = 1) - (Z[i] == 0) / GPS(z = 0)

    if (length(neighbor) > 0) {

      W_Y[i]      <- sum(Y0[neighbor]) * multiplier
      W_D_AIED[i] <- sum(D0[neighbor]) * multiplier
      W_D_ADED[i] <- D[i] * multiplier

    }

  }

  V_AIEY <- W_Y -
    W_Z1 * mu_Y(z = 1) / GPS(z = 1) +
    W_Z0 * mu_Y(z = 0) / GPS(z = 0)

  V_AIED <- W_D_AIED -
    W_Z1 * mu_D_AIED(z = 1) / GPS(z = 1) +
    W_Z0 * mu_D_AIED(z = 0) / GPS(z = 0)

  V_ADED <- W_D_ADED -
    W_Z1 * mu_D_ADED(z = 1) / GPS(z = 1) +
    W_Z0 * mu_D_ADED(z = 0) / GPS(z = 0)

  V_LAIE <- V_AIEY / est_ADED - V_ADED * est_AIEY / est_ADED^2

  # Function to compute HAC standard error and confidence interval
  # Input: Nothing
  # Output: A list containing HAC standard error and confidence interval
  HAC_func <- function() {

    # Standard error
    neighbor_mat <- (distance <= bw)

    SE_AIEY <- sqrt(t(V_AIEY) %*% neighbor_mat %*% V_AIEY / size^2)
    SE_AIED <- sqrt(t(V_AIED) %*% neighbor_mat %*% V_AIED / size^2)
    SE_ADED <- sqrt(t(V_ADED) %*% neighbor_mat %*% V_ADED / size^2)
    SE_LAIE <- sqrt(t(V_LAIE) %*% neighbor_mat %*% V_LAIE / size^2)

    # Critical value based on asymptotic normality
    cvnorm <- stats::qnorm(1 - alp / 2)

    # Confidence interval
    CI_AIEY <- c(est_AIEY - cvnorm * SE_AIEY,
                 est_AIEY + cvnorm * SE_AIEY)

    CI_AIED <- c(est_AIED - cvnorm * SE_AIED,
                 est_AIED + cvnorm * SE_AIED)

    CI_ADED <- c(est_ADED - cvnorm * SE_ADED,
                 est_ADED + cvnorm * SE_ADED)

    CI_LAIE <- c(est_LAIE - cvnorm * SE_LAIE,
                 est_LAIE + cvnorm * SE_LAIE)

    # Result
    HAC_AIEY <- c(SE_AIEY, CI_AIEY)
    HAC_AIED <- c(SE_AIED, CI_AIED)
    HAC_ADED <- c(SE_ADED, CI_ADED)
    HAC_LAIE <- c(SE_LAIE, CI_LAIE)

    return(list(AIEY = HAC_AIEY,
                AIED = HAC_AIED,
                ADED = HAC_ADED,
                LAIE = HAC_LAIE))

  }

  # Function to compute the Omega matrix used for wild bootstrap
  # Input: Nothing
  # Output: the Omega matrix
  Omega_func <- function() {

    # Compute the numerator
    neighbor_mat <- (distance <= bw)
    numerator <- NULL
    for (i in 1:size) {

      intersect_mat <- sweep(x = neighbor_mat,
                             MARGIN = 2,
                             FUN = "*",
                             neighbor_mat[i, ])

      numerator <- rbind(numerator,
                         apply(intersect_mat, MARGIN = 1, FUN = sum))

    }

    # Compute the denominator
    denominator <- sum(neighbor_mat) / size

    return(numerator/denominator)

  }

  # Function to compute wild bootstrap standard error and confidence interval
  # Input: Nothing
  # Output: A list containing bootstrap standard error and confidence interval
  wild_func <- function() {

    # square root matrix of Omega
    Omega <- Omega_func()
    Omega <- methods::as(Omega, "dgCMatrix")
    E <- eigen(Omega)
    E$values[abs(E$values) < 1e-10] <- 0
    sqrt_Omega <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)

    # bootstrap
    AIEY_boot <- AIED_boot <- ADED_boot <- LAIE_boot <- NULL
    for (r in 1:B) {

      zeta <- stats::rnorm(size)

      R_vector <- sqrt_Omega %*% zeta

      V_AIEY_boot <- V_AIEY * R_vector
      V_AIED_boot <- V_AIED * R_vector
      V_ADED_boot <- V_ADED * R_vector
      V_LAIE_boot <- V_LAIE * R_vector

      AIEY_boot <- c(AIEY_boot, mean(V_AIEY_boot))
      AIED_boot <- c(AIED_boot, mean(V_AIED_boot))
      ADED_boot <- c(ADED_boot, mean(V_ADED_boot))
      LAIE_boot <- c(LAIE_boot, mean(V_LAIE_boot))

    }

    # standard error
    SE_AIEY <- stats::sd(AIEY_boot)
    SE_AIED <- stats::sd(AIED_boot)
    SE_ADED <- stats::sd(ADED_boot)
    SE_LAIE <- stats::sd(LAIE_boot)

    # confidence interval
    q_AIEY <- stats::quantile(AIEY_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    q_AIED <- stats::quantile(AIED_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    q_ADED <- stats::quantile(ADED_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    q_LAIE <- stats::quantile(LAIE_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    CI_AIEY <- c(est_AIEY + q_AIEY[1],
                 est_AIEY + q_AIEY[2])

    CI_AIED <- c(est_AIED + q_AIED[1],
                 est_AIED + q_AIED[2])

    CI_ADED <- c(est_ADED + q_ADED[1],
                 est_ADED + q_ADED[2])

    CI_LAIE <- c(est_LAIE + q_LAIE[1],
                 est_LAIE + q_LAIE[2])

    # Result
    wild_AIEY <- c(SE_AIEY, CI_AIEY)
    wild_AIED <- c(SE_AIED, CI_AIED)
    wild_ADED <- c(SE_ADED, CI_ADED)
    wild_LAIE <- c(SE_LAIE, CI_LAIE)

    return(list(AIEY = wild_AIEY,
                AIED = wild_AIED,
                ADED = wild_ADED,
                LAIE = wild_LAIE))

  }

  # HAC estimation
  HAC <- HAC_func()

  # Wild bootstrap
  if (!is.null(B)) {

    wild <- wild_func()

  } else {

    wild <- list(AIEY = rep(NA, 3),
                 AIED = rep(NA, 3),
                 ADED = rep(NA, 3),
                 LAIE = rep(NA, 3))

  }

  # Result
  Estimate <- rbind(c(est_AIEY, HAC$AIEY, wild$AIEY, bw, size),
                    c(est_AIED, HAC$AIED, wild$AIED, bw, size),
                    c(est_ADED, HAC$ADED, wild$ADED, bw, size),
                    c(est_LAIE, HAC$LAIE, wild$LAIE, bw, size)
  )

  colnames(Estimate) <- c("est",
                          "HAC_SE",
                          "HAC_CI_L",
                          "HAC_CI_U",
                          "wild_SE",
                          "wild_CI_L",
                          "wild_CI_U",
                          "bw",
                          "size")

  rownames(Estimate) <- c("AIEY", "AIED", "ADED", "LAIE")

  Estimate <- as.data.frame(Estimate)

  return(Estimate)

}
