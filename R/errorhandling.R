#' Error handling
#'
#' @param Y An n-dimensional outcome vector
#' @param D An n-dimensional binary treatment vector
#' @param Z An n-dimensional binary instrumental vector
#' @param IEM An n-dimensional instrumental exposure vector
#' @param S An n-dimensional logical vector of indicating sub-population S
#' @param A An n times n symmetric binary adjacency matrix
#' @param bw A scalar of bandwidth used for HAC estimation and wild bootstrap
#' @param B The number of bootstrap repetitions
#' @param alp The significance level
#'
#' @returns NULL if there are no errors
#'
#' @noRd
#'
errorhandling <- function(Y,
                          D,
                          Z,
                          IEM,
                          S,
                          A,
                          bw,
                          B,
                          alp) {

  # Y --------------------------------------------------------------------------

  if (!is.numeric(Y)) {

    # Numerical vector
    stop(paste("Y must be a numerical vector."))

  } else {

    # NA
    if (any(is.na(Y))) {
      stop(paste("Y must not contain NA."))
    }

  }

  # D --------------------------------------------------------------------------

  if (!is.numeric(D)) {

    # Numerical vector
    stop(paste("D must be a numerical vector."))

  } else {

    # NA
    if (any(is.na(D))) {
      stop(paste("D must not contain NA."))
    }

    # Binary
    if (sum(D != 0 & D != 1) > 0) {
      stop(paste("D must not contain the values other than 0 and 1."))
    }

  }

  # Z --------------------------------------------------------------------------

  if (!is.numeric(Z)) {

    # Numerical vector
    stop(paste("Z must be a numerical vector."))

  } else {

    # NA
    if (any(is.na(Z))) {
      stop(paste("Z must not contain NA."))
    }

    # Binary
    if (sum(Z != 0 & Z != 1) > 0) {
      stop(paste("Z must not contain the values other than 0 and 1."))
    }

  }

  # IEM ------------------------------------------------------------------------

  if (!is.null(IEM)) {

    if (!is.numeric(IEM)) {

      # Numerical vector
      stop(paste("IEM must be a numerical vector or NULL."))

    } else {

      # NA
      if (any(is.na(IEM))) {
        stop(paste("IEM must not contain NA."))
      }

    }

  }

  # S --------------------------------------------------------------------------

  if (!is.logical(S)) {

    # Logical vector
    stop(paste("S must be a logical vector."))

  } else {

    # NA
    if (any(is.na(S))) {
      stop(paste("S must not contain NA."))
    }

  }

  # A --------------------------------------------------------------------------

  if (!is.matrix(A)) {

    # Matrix
    stop(paste("A must be a symmetric binary matrix."))

  } else {

    # NA
    if (any(is.na(A))) {
      stop(paste("A must not contain NA."))
    }

    # Binary matrix
    if (sum(A != 0 & A != 1) > 0) {
      stop(paste("A must not contain the values other than 0 and 1."))
    }

    # Symmetric matrix
    if (sum(A != t(A)) != 0) {
      stop(paste("A must be a symmetric binary matrix."))
    }

    # Diagonal elements
    if (sum(diag(A) != 0) > 0) {
      stop(paste("The diagonal elements of A must be 0."))
    }

    # Square matrix
    if (nrow(A) != ncol(A)) {
      stop(paste("A must be a square matrix."))
    }

  }

  # Dimensions -----------------------------------------------------------------

  dim_Y <- length(Y)
  dim_D <- length(D)
  dim_Z <- length(Z)
  dim_S <- length(S)
  dim_A <- nrow(A)

  if (dim_Y != dim_D | dim_Y != dim_Z | dim_Y != dim_S | dim_Y != dim_A) {
    stop(paste("The lengths of Y, D, Z, S,
               and of the row and column of A must be the same."))
  }

  # bw -------------------------------------------------------------------------

  if (!is.null(bw)) {

    if (!is.numeric(bw)) {

      # Numeric
      stop(paste("bw must be a non-negative number."))

    } else {

      # Scalar
      if (length(bw) != 1) {
        stop(paste("bw must be a non-negative number."))
      }

      # Non-negative number
      if (bw < 0) {
        stop(paste("bw must be a non-negative number."))
      }

    }

  }

  # B --------------------------------------------------------------------------

  if (!is.null(B)) {

    # Numeric
    if (!is.numeric(B)) {

      stop(paste("B must be a positive number."))

    } else {

      # Scalar
      if (length(B) != 1) {
        stop(paste("B must be a positive number."))
      }

      # Positive number
      if (B <= 0) {
        stop(paste("B must be a positive number."))
      }

    }

  }

  # alp ------------------------------------------------------------------------

  if (!is.numeric(alp)) {

    # Numeric
    stop(paste("alp must be a positive number between 0 and 1."))

  } else {

    # Scalar
    if (length(alp) != 1) {
      stop(paste("alp must be a positive number between 0 and 1."))
    }

    # Positive value between 0 and 1
    if (alp <= 0 | alp >= 1) {
      stop(paste("alp must be a positive number between 0 and 1."))
    }

  }

  # Return ---------------------------------------------------------------------

  return(NULL)

}
