
<!-- README.md is generated from README.Rmd. Please edit that file -->

# latenetwork: Inference on LATEs under Network Interference of Unknown Form

<!-- badges: start -->
<!-- badges: end -->

The **latenetwork** package provides tools for causal inference under
noncompliance with treatment assignment and network interference of
unknown form. The package enables to implement the instrumental
variables (IV) estimation for the local average treatment effect (LATE)
type parameters via inverse probability weighting (IPW) using the
concept of instrumental exposure mapping (IEM) and the framework of
approximate neighborhood interference (ANI). For more details, see
[Hoshino and Yanagi (2023) “Causal inference with noncompliance and
unknown interference”](https://doi.org/10.48550/arXiv.2108.07455).

## Installation

Get the package from CRAN:

``` r
install.packages("latenetwork")
```

or from GitHub:

``` r
# install.packages("devtools") # if needed
devtools::install_github("tkhdyanagi/latenetwork", build_vignettes = TRUE)
```

## Vignettes

For more details, see the package vignettes with:

``` r
library("latenetwork")

# Getting Started with the latenetwork Package
vignette("latenetwork")

# Review of Causal Inference with Noncompliance and Unknown Interference
vignette("review", package = "latenetwork")
```

## References

- Hoshino, T. and Yanagi, T., 2023. Causal inference with noncompliance
  and unknown interference. arXiv preprint arXiv:2108.07455.
  [Link](https://doi.org/10.48550/arXiv.2108.07455)
