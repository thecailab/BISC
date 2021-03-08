#' The 'BISC' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
#' Gene expression in mammalian cells is inherently stochastic that mRNA molecules are synthesized in discrete bursts. Although the advent of single cell transcriptomic technology provides great opportunities to explore the phenomenon of transcriptional bursting, current Beta-Poisson framework greatly suffers from substantial technical noise that leads to false estimations and conclusions. Statistical methods that account for complex nature of single cell transcriptomic data are needed to accurately reveal the heterogeneity of gene expression and infer kinetics parameters. Here, we develop a Bayesian hierarchical framework, BISC, to study the stochastic gene expression kinetics. Gamma-Poisson model was employed to accommodate overdispersion of read counts, which was dynamically characterized by fitting the mean-variance relationship. The reliable estimation of the dispersion parameter is essential for capturing the expression heterogeneity, which improve the estimation of kinetics parameter. Also, we proposed a differential bursting analysis framework to identify heterogeneous bursting kinetics under different studying conditions.
#' @docType package
#' @name BISC-package
#' @aliases BISC
#' @useDynLib BISC, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.1. https://mc-stan.org
#'
NULL
