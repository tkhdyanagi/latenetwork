#' Inference on Average Direct Effect Parameters
#'
#' Inference on the average direct effect of the IV on the outcome,
#' that on the treatment receipt, and the local average direct effect
#' in the presence of network spillover of unknown form
#'
#' @details
#' The `direct()` function estimates the average direct effect of the IV
#' on the outcome, that on the treatment receipt, and
#' the local average direct effect via inverse probability weighting
#' in the approximate neighborhood interference framework.
#' The function also computes the standard errors and the confidence intervals
#' for the target parameters based on the network HAC estimation and
#' the wild bootstrap.
#' For more details, see Hoshino and Yanagi (2023).
#' The lengths of `Y`, `D`, `Z`, `S` and
#' of the row and column of `A` must be the same.
#' `IEM` must be `NULL` or a vector of the same length as `Y`.
#' `t` must be `NULL` or a value in the support of `IEM`.
#' `K` must be a positive integer.
#' `bw` must be `NULL` or a non-negative integer.
#' `B` must be `NULL` or a positive number.
#' `alp` must be a positive number between 0 and 0.5.
#'
#' @param Y An n-dimensional outcome vector
#' @param D An n-dimensional binary treatment vector
#' @param Z An n-dimensional binary instrumental vector
#' @param IEM An n-dimensional instrumental exposure vector.
#' If `IEM = NULL` or `t = NULL`, the constant IEM is used.
#' Default is NULL.
#' @param S An n-dimensional logical vector to indicate whether each unit
#' belongs to the sub-population S
#' @param A An n times n symmetric binary adjacency matrix
#' @param K A scalar to indicate the range of neighborhood
#' used for constructing the interference set.
#' Default is 1.
#' In the `direct()` function, `K` is used only for computing the bandwidth.
#' @param t A scalar of the evaluation point of IEM.
#' Default is NULL.
#' @param bw A scalar of the bandwidth used for the HAC estimation and
#' the wild bootstrap.
#' If `bw = NULL`, the rule-of-thumb bandwidth proposed by Leung (2022) is used.
#' Default is NULL.
#' @param B The number of bootstrap repetitions.
#' If `B = NULL`, the wild bootstrap is skipped.
#' Default is NULL.
#' @param alp The significance level.
#' Default is 0.05.
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
#' IEM <- data$IEM
#' S   <- rep(TRUE, n)
#' A   <- data$A
#' K   <- 1
#' t   <- 0
#' bw  <- NULL
#' B   <- NULL
#' alp <- 0.05
#'
#' # Estimation
#' latenetwork::direct(Y = Y,
#'                     D = D,
#'                     Z = Z,
#'                     IEM = IEM,
#'                     S = S,
#'                     A = A,
#'                     K = K,
#'                     t = t,
#'                     bw = bw,
#'                     B = B,
#'                     alp = alp)
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
direct <- function(Y,
                   D,
                   Z,
                   IEM = NULL,
                   S,
                   A,
                   K = 1,
                   t = NULL,
                   bw = NULL,
                   B = NULL,
                   alp = 0.05) {

  # Error handling -------------------------------------------------------------

  error <- errorhandling(Y = Y,
                         D = D,
                         Z = Z,
                         IEM = IEM,
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
                                 v =  igraph::V(graph0),
                                 to = igraph::V(graph0),
                                 algorithm = "dijkstra")

  # Variable definitions on the sub-population S
  Y <- Y[S]
  D <- D[S]
  Z <- Z[S]
  A <- A[S, S]
  graph <- igraph::graph_from_adjacency_matrix(adjmatrix = A,
                                               mode = "undirected")
  distance <- distance0[S, S]

  # Constant IEM or not
  if (is.null(IEM) | is.null(t)) {
    IEM <- t <- TRUE
  } else {
    IEM <- IEM[S]
  }

  # Size of the sub-population S
  size <- sum(S)

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
  # Input "t": A value of IEM
  # Output: A value of GPS
  GPS <- function(z, t) {
    mean(Z == z & IEM == t)
  }

  # Function to compute the conditional outcome mean
  # Input "z": A value of IV
  # Input "t": A value of IEM
  # Output: A value of the conditional outcome mean
  mu_Y <- function(z, t) {
    mean(Y * (Z == z & IEM == t)) / GPS(z = z, t = t)
  }

  # Function to compute the conditional treatment mean
  # Input "z": A value of IV
  # Input "t": A value of IEM
  # Output: A value of the conditional treatment mean
  mu_D <- function(z, t) {
    mean(D * (Z == z & IEM == t)) / GPS(z = z, t = t)
  }

  # Function to compute ADEY(t)
  # Input "t": A value of IEM
  # Output: A value of ADEY(t)
  ADEY <- function(t) {
    mu_Y(z = 1, t = t) - mu_Y(z = 0, t = t)
  }

  # Function to compute ADED(t)
  # Input "t": A value of IEM
  # Output: A value of ADED(t)
  ADED <- function(t) {
    mu_D(z = 1, t = t) - mu_D(z = 0, t = t)
  }

  # Function to compute LADE(t)
  # Input "t": A value of IEM
  # Output: A value of LADE(t)
  LADE <- function(t) {
    ADEY(t = t) / ADED(t = t)
  }

  # Estimates
  est_ADEY <- ADEY(t = t)
  est_ADED <- ADED(t = t)
  est_LADE <- LADE(t = t)

  # Variable definitions
  W_Z1 <- (Z == 1 & IEM == t)

  W_Z0 <- (Z == 0 & IEM == t)

  W_Y <- Y * (W_Z1 / GPS(z = 1, t = t) - W_Z0 / GPS(z = 0, t = t))

  W_D <- D * (W_Z1 / GPS(z = 1, t = t) - W_Z0 / GPS(z = 0, t = t))

  V_ADEY <- W_Y -
    W_Z1 * mu_Y(z = 1, t = t) / GPS(z = 1, t = t) +
    W_Z0 * mu_Y(z = 0, t = t) / GPS(z = 0, t = t)

  V_ADED <- W_D -
    W_Z1 * mu_D(z = 1, t = t) / GPS(z = 1, t = t) +
    W_Z0 * mu_D(z = 0, t = t) / GPS(z = 0, t = t)

  V_LADE <- V_ADEY / est_ADED - V_ADED * est_ADEY / est_ADED^2

  # Function to compute HAC standard error and confidence interval
  # Input: Nothing
  # Output: A list containing HAC standard error and confidence interval
  HAC_func <- function() {

    # HAC standard error
    neighbor_mat <- (distance <= bw)

    SE_ADEY <- sqrt(t(V_ADEY) %*% neighbor_mat %*% V_ADEY / size^2)
    SE_ADED <- sqrt(t(V_ADED) %*% neighbor_mat %*% V_ADED / size^2)
    SE_LADE <- sqrt(t(V_LADE) %*% neighbor_mat %*% V_LADE / size^2)

    # Critical value based on asymptotic normality
    cvnorm <- stats::qnorm(1 - alp / 2)

    # Confidence interval
    CI_ADEY <- c(est_ADEY - cvnorm * SE_ADEY,
                 est_ADEY + cvnorm * SE_ADEY)

    CI_ADED <- c(est_ADED - cvnorm * SE_ADED,
                 est_ADED + cvnorm * SE_ADED)

    CI_LADE <- c(est_LADE - cvnorm * SE_LADE,
                 est_LADE + cvnorm * SE_LADE)

    # results
    HAC_ADEY <- c(SE_ADEY, CI_ADEY)
    HAC_ADED <- c(SE_ADED, CI_ADED)
    HAC_LADE <- c(SE_LADE, CI_LADE)

    return(list(ADEY = HAC_ADEY,
                ADED = HAC_ADED,
                LADE = HAC_LADE))

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

    # Square root matrix of Omega
    Omega <- Omega_func()
    Omega <- methods::as(Omega, "dgCMatrix")
    E <- eigen(Omega)
    E$values[abs(E$values) < 1e-10] <- 0
    sqrt_Omega <- E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)

    # Bootstrap
    ADED_boot <- ADEY_boot <- LADE_boot <- NULL
    for (r in 1:B) {

      zeta <- stats::rnorm(size)

      R_vector <- sqrt_Omega %*% zeta

      V_ADEY_boot <- V_ADEY * R_vector
      V_ADED_boot <- V_ADED * R_vector
      V_LADE_boot <- V_LADE * R_vector

      ADEY_boot <- c(ADEY_boot, mean(V_ADEY_boot))
      ADED_boot <- c(ADED_boot, mean(V_ADED_boot))
      LADE_boot <- c(LADE_boot, mean(V_LADE_boot))

    }

    # Standard error
    SE_ADEY <- stats::sd(ADEY_boot)
    SE_ADED <- stats::sd(ADED_boot)
    SE_LADE <- stats::sd(LADE_boot)

    # Confidence interval
    q_ADEY <- stats::quantile(ADEY_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    q_ADED <- stats::quantile(ADED_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    q_LADE <- stats::quantile(LADE_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    CI_ADEY <- c(est_ADEY + q_ADEY[1],
                 est_ADEY + q_ADEY[2])

    CI_ADED <- c(est_ADED + q_ADED[1],
                 est_ADED + q_ADED[2])

    CI_LADE <- c(est_LADE + q_LADE[1],
                 est_LADE + q_LADE[2])

    # Result
    wild_ADEY <- c(SE_ADEY, CI_ADEY)
    wild_ADED <- c(SE_ADED, CI_ADED)
    wild_LADE <- c(SE_LADE, CI_LADE)

    return(list(ADEY = wild_ADEY,
                ADED = wild_ADED,
                LADE = wild_LADE))

  }

  # HAC estimation
  HAC <- HAC_func()

  # Wild bootstrap
  if (!is.null(B)) {

    wild <- wild_func()

  } else {

    wild <- list(ADEY = rep(NA, 3),
                 ADED = rep(NA, 3),
                 LADE = rep(NA, 3))

  }

  # Results
  Estimate <- rbind(c(est_ADEY, HAC$ADEY, wild$ADEY, bw, size),
                    c(est_ADED, HAC$ADED, wild$ADED, bw, size),
                    c(est_LADE, HAC$LADE, wild$LADE, bw, size)
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

  rownames(Estimate) <- c("ADEY", "ADED", "LADE")

  Estimate <- as.data.frame(Estimate)

  return(Estimate)

}
