#' Inference on Average Spillover Effect Parameters
#'
#' Inference on the average spillover effect of the IV on the outcome,
#' that on the treatment receipt, and the local average spillover effect
#' in the presence of network spillover of unknown form
#'
#' @details
#' The `spillover()` function estimates the average spillover effect of the IV
#' on the outcome, that on the treatment receipt,
#' and the local average spillover effect via inverse probability weighting
#' in the approximate neighborhood interference framework.
#' The function also computes the standard errors and the confidence intervals
#' for the target parameters based on the network HAC estimation and
#' the wild bootstrap.
#' For more details, see Hoshino and Yanagi (2023).
#' The lengths of `Y`, `D`, `Z`, `IEM`, `S` and
#' of the row and column of `A` must be the same.
#' `z` must be 0 or 1.
#' `t0` and `t1` must be values in the support of `IEM`.
#' `bw` must be `NULL` or a non-negative number.
#' `B` must be `NULL` or a positive integer.
#' `alp` must be a positive number between 0 and 0.5.
#'
#' @param Y An n-dimensional outcome vector
#' @param D An n-dimensional binary treatment vector
#' @param Z An n-dimensional binary instrumental vector
#' @param IEM An n-dimensional instrumental exposure vector
#' @param S An n-dimensional logical vector to indicate whether each unit
#' belongs to the sub-population S
#' @param A An n times n symmetric binary adjacency matrix
#' @param K A scalar to indicate the range of neighborhood
#' used for constructing the interference set.
#' Default is 1.
#' In the `spillover()` function, `K` is used only for computing the bandwidth.
#' @param z A scalar of the evaluation point of Z
#' @param t0 A scalar of the evaluation point of instrumental exposure (from)
#' @param t1 A scalar of the evaluation point of instrumental exposure (to)
#' @param bw A scalar of the bandwidth used for the HAC estimation and
#' the wild bootstrap.
#' If `bw = NULL`, the rule-of-thumb bandwidth proposed by Leung (2022) is used.
#' Default is NULL.
#' @param B The number of bootstrap repetitions.
#' If `B = NULL`, wild bootstrap is skipped.
#' Default is NULL.
#' @param alp The significance level.
#' Default is 0.05.
#'
#' @returns A data frame containing the following elements:
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
#' IEM <- ifelse(A %*% Z > 0, 1, 0)
#' z   <- 1
#' t0  <- 0
#' t1  <- 1
#' bw  <- NULL
#' B   <- NULL
#' alp <- 0.05
#'
#' # Estimation
#' latenetwork::spillover(Y = Y,
#'                        D = D,
#'                        Z = Z,
#'                        IEM = IEM,
#'                        S = S,
#'                        A = A,
#'                        K = K,
#'                        z = z,
#'                        t0 = t0,
#'                        t1 = t1,
#'                        bw = bw,
#'                        B = B,
#'                        alp = alp)
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
spillover <- function(Y,
                      D,
                      Z,
                      IEM,
                      S,
                      A,
                      K = 1,
                      z,
                      t0,
                      t1,
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
                                 v  = igraph::V(graph0),
                                 to = igraph::V(graph0),
                                 algorithm = "dijkstra")

  # Variable definitions on the sub-population S
  Y <- Y[S]
  D <- D[S]
  Z <- Z[S]
  IEM <- IEM[S]
  A <- A[S, S]
  graph <- igraph::graph_from_adjacency_matrix(adjmatrix = A,
                                               mode = "undirected")
  distance <- distance0[S, S]

  # Sample size
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

  # Function to compute ASEY(z, t0, t1)
  # Input "z": A value of IV
  # Input "t0": A value of IEM (from)
  # Input "t1": A value of IEM (to)
  # Output: A value of ASEY(z, t0, t1)
  ASEY <- function(z, t0, t1) {
    mu_Y(z = z, t = t1) - mu_Y(z = z, t = t0)
  }

  # Function to compute ASED(z, t0, t1)
  # Input "z": A value of IV
  # Input "t0": A value of IEM (from)
  # Input "t1": A value of IEM (to)
  # Output: A value of ASED(z, t0, t1)
  ASED <- function(z, t0, t1) {
    mu_D(z = z, t = t1) - mu_D(z = z, t = t0)
  }

  # Function to compute LASE(z, t0, t1)
  # Input "z": A value of IV
  # Input "t0": A value of IEM (from)
  # Input "t1": A value of IEM (to)
  # Output: A value of LASE(z, t0, t1)
  LASE <- function(z, t0, t1) {
    ASEY(z = z, t0 = t0, t1 = t1) / ASED(z = z, t0 = t0, t1 = t1)
  }

  # Estimates
  est_ASEY <- ASEY(z = z, t0 = t0, t1 = t1)
  est_ASED <- ASED(z = z, t0 = t0, t1 = t1)
  est_LASE <- LASE(z = z, t0 = t0, t1 = t1)

  # Variable definitions
  X_IEM1 <- (Z == z & IEM == t1)

  X_IEM0 <- (Z == z & IEM == t0)

  X_Y <- Y * (X_IEM1 / GPS(z = z, t = t1) - X_IEM0 / GPS(z = z, t = t0))

  X_D <- D * (X_IEM1 / GPS(z = z, t = t1) - X_IEM0 / GPS(z = z, t = t0))

  V_ASEY <- X_Y -
    X_IEM1 * mu_Y(z = z, t = t1) / GPS(z = z, t = t1) +
    X_IEM0 * mu_Y(z = z, t = t0) / GPS(z = z, t = t0)

  V_ASED <- X_D -
    X_IEM1 * mu_D(z = z, t = t1) / GPS(z = z, t = t1) +
    X_IEM0 * mu_D(z = z, t = t0) / GPS(z = z, t = t0)

  V_LASE <- V_ASEY / est_ASED - V_ASED * est_ASEY / est_ASED^2

  # Function to compute HAC standard error and confidence interval
  # Input: Nothing
  # Output: A list containing HAC standard error and confidence interval
  HAC_func <- function() {

    # HAC standard error
    neighbor_mat <- (distance <= bw)

    SE_ASEY <- sqrt(t(V_ASEY) %*% neighbor_mat %*% V_ASEY / size^2)
    SE_ASED <- sqrt(t(V_ASED) %*% neighbor_mat %*% V_ASED / size^2)
    SE_LASE <- sqrt(t(V_LASE) %*% neighbor_mat %*% V_LASE / size^2)

    # Critical value based on asymptotic normality
    cvnorm <- stats::qnorm(1 - alp / 2)

    # Confidence interval
    CI_ASEY <- c(est_ASEY - cvnorm * SE_ASEY,
                 est_ASEY + cvnorm * SE_ASEY)

    CI_ASED <- c(est_ASED - cvnorm * SE_ASED,
                 est_ASED + cvnorm * SE_ASED)

    CI_LASE <- c(est_LASE - cvnorm * SE_LASE,
                 est_LASE + cvnorm * SE_LASE)

    # Results
    HAC_ASEY <- c(SE_ASEY, CI_ASEY)
    HAC_ASED <- c(SE_ASED, CI_ASED)
    HAC_LASE <- c(SE_LASE, CI_LASE)

    return(list(ASEY = HAC_ASEY,
                ASED = HAC_ASED,
                LASE = HAC_LASE))

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
    ASEY_boot <- ASED_boot <- LASE_boot <- NULL
    for (r in 1:B) {

      zeta <- stats::rnorm(size)

      R_vector <- sqrt_Omega %*% zeta

      V_ASEY_boot <- V_ASEY * R_vector
      V_ASED_boot <- V_ASED * R_vector
      V_LASE_boot <- V_LASE * R_vector

      ASEY_boot <- c(ASEY_boot, mean(V_ASEY_boot))
      ASED_boot <- c(ASED_boot, mean(V_ASED_boot))
      LASE_boot <- c(LASE_boot, mean(V_LASE_boot))

    }

    # Standard error
    SE_ASEY <- stats::sd(ASEY_boot)
    SE_ASED <- stats::sd(ASED_boot)
    SE_LASE <- stats::sd(LASE_boot)

    # Confidence interval
    q_ASEY <- stats::quantile(ASEY_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    q_ASED <- stats::quantile(ASED_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    q_LASE <- stats::quantile(LASE_boot,
                              probs = c(alp / 2, 1 - alp / 2))

    CI_ASEY <- c(est_ASEY + q_ASEY[1],
                 est_ASEY + q_ASEY[2])

    CI_ASED <- c(est_ASED + q_ASED[1],
                 est_ASED + q_ASED[2])

    CI_LASE <- c(est_LASE + q_LASE[1],
                 est_LASE + q_LASE[2])

    # Results
    wild_ASEY <- c(SE_ASEY, CI_ASEY)
    wild_ASED <- c(SE_ASED, CI_ASED)
    wild_LASE <- c(SE_LASE, CI_LASE)

    return(list(ASEY = wild_ASEY,
                ASED = wild_ASED,
                LASE = wild_LASE))

  }

  # HAC estimation
  HAC <- HAC_func()

  # Wild bootstrap
  if (!is.null(B)) {

    wild <- wild_func()

  } else {

    wild <- list(ASEY = rep(NA, 3),
                 ASED = rep(NA, 3),
                 LASE = rep(NA, 3))

  }

  # Results
  Estimate <- rbind(
    c(est_ASEY, HAC$ASEY, wild$ASEY, bw, size),
    c(est_ASED, HAC$ASED, wild$ASED, bw, size),
    c(est_LASE, HAC$LASE, wild$LASE, bw, size)
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

  rownames(Estimate) <- c("ASEY", "ASED", "LASE")

  Estimate <- as.data.frame(Estimate)

  return(Estimate)

}
