# Simulated data ---------------------------------------------------------------
set.seed(1)

# Sample size
n <- 1000

# Adjacency matrix of ring network
ring <- igraph::make_ring(n, directed = FALSE, mutual = FALSE, circular = TRUE)
A <- igraph::as_adjacency_matrix(ring)
A <- as.matrix(A)

# Instrumental vector
Z <- sample(x = 0:1, size = n, replace = TRUE, prob = c(0.5, 0.5))

# Instrumental exposure
IEM <- ifelse(A %*% Z > 0, 1, 0)

# Treatment vector
D <- ifelse(rnorm(n, mean = -1.5, sd = 1) + Z + IEM >= 0, 1, 0)

# Outcome vector
Y <- rnorm(n, mean = 1, sd = 1) + runif(n, min = 1, max = 2) * D

# Other arguments
S <- rep(TRUE, n)
K <- 1
z <- 1
t <- 0
t0 <- 0
t1 <- 1
bw <- NULL
B <- NULL
alp <- 0.05

# Y ----------------------------------------------------------------------------

expect_error(direct(Y = TRUE,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "Y must be a numerical vector.")

Y1 <- Y
Y1[1] <- NA

expect_error(direct(Y = Y1,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "Y must not contain NA.")

# D ----------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = FALSE,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "D must be a numerical vector.")

D1 <- D
D1[2] <- NA

expect_error(direct(Y = Y,
                    D = D1,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "D must not contain NA.")

D1 <- D
D1[3] <- 5

expect_error(direct(Y = Y,
                    D = D1,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "D must not contain the values other than 0 and 1.")

# Z ----------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = D,
                    Z = NULL,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "Z must be a numerical vector.")

Z1 <- Z
Z1[3] <- NA

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z1,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "Z must not contain NA.")

Z1 <- Z
Z1[4] <- 100

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z1,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "Z must not contain the values other than 0 and 1.")

# IEM --------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = Sys.time(),
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "IEM must be a numerical vector or NULL.")

IEM1 <- IEM
IEM1[10] <- NA

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM1,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "IEM must not contain NA.")

# S ----------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = 1:n,
                    A = A,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "S must be a logical vector.")

S1 <- S
S1[100] <- NA

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S1,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "S must not contain NA.")

# A ----------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = NULL,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "A must be a symmetric binary matrix.")

A1 <- A
A1[1, 1] <- NA

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A1,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "A must not contain NA.")

A2 <- A
A2[1, 20] <- 7

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A2,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "A must not contain the values other than 0 and 1.")

A3 <- A
A3[1, 2] <- 0

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A3,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "A must be a symmetric binary matrix.")

A4 <- A
A4[1, 1] <- 1

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A4,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "The diagonal elements of A must be 0.")

A5 <- A[-1, ]

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A5,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp))

# Dimensions -------------------------------------------------------------------

Y2 <- Y[-1]

expect_error(direct(Y = Y2,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = alp),
             "The lengths of Y, D, Z, S,
               and of the row and column of A must be the same.")

# bw ---------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = NA,
                    B = B,
                    alp = alp),
             "bw must be a non-negative number.")

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = 1:2,
                    B = B,
                    alp = alp),
             "bw must be a non-negative number.")

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = -1,
                    B = B,
                    alp = alp),
             "bw must be a non-negative number.")

# B ----------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = TRUE,
                    alp = alp),
             "B must be a positive number.")

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = 1:n,
                    alp = alp),
             "B must be a positive number.")

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = -100,
                    alp = alp),
             "B must be a positive number.")

# alp --------------------------------------------------------------------------

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = NULL),
             "alp must be a positive number between 0 and 1.")

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = c(0.5, 0.95)),
             "alp must be a positive number between 0 and 1.")

expect_error(direct(Y = Y,
                    D = D,
                    Z = Z,
                    IEM = IEM,
                    S = S,
                    A = A,
                    K = K,
                    t = t,
                    bw = bw,
                    B = B,
                    alp = 1.5),
             "alp must be a positive number between 0 and 1.")
