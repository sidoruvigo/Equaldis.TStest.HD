#' Package ‘Equaldis.TStest.HD’
#'
#' Documentation for package ‘Equaldis.TStest.HD’ version 1.0
#' @aliases Equaldis.TStest.HDpackage
#' @title A two-sample test for the equality of distributions for high-dimensional data
#'
#' @description This package implements four different methods to test the (global) null hypothesis of equality of the
#'  p marginal distributions for two groups which have been proposed in Cousido-Rocha et al. (2018).
#'  In other words, the null hypothesis is an intersection of the p null hypotheses corresponding to p different
#'  two-sample problems. This null hypothesis implies that the univariate marginals of the p-variate distributions
#'  coincide in the two populations. These methods are particularly well suited to the low sample size, high
#'  dimensional setting (n << p). The sample size can be as small as 2. The test accounts for the possibility
#'  that the p variables in each data set can be weakly dependent. Three of the methods arise from different
#'  approaches to estimate the variance of the same statistic. This statistic is based on the simple idea of
#'  comparing, for each of the p variables, the empirical characteristic functions computed from the two
#'  samples. The last method is an alternative global test based on the p-values derived from permutation
#'  tests. When the global null hypothesis is rejected such permutation p-values can be used to identify
#'  which variables have been contributed to this significance. The standarized version of each test statistic
#'  and its p-value are computed among other things.
#'
#'
#' @details
#' \itemize{
#' \item{Package: ‘Equaldis.TStest.HD’}
#' \item{Version: 1.0}
#' \item{Maintainer: Marta Cousido Rocha \email{martacousido@@uvigo.es}}
#' \item{License: GPL-2}
#' }
#'
#' @return
#' \itemize{
#' \item{‘Equaldis.TStest.HD’}
#' }
#'
#' @author
#' \itemize{
#' \item{Cousido Rocha, Marta.}
#' \item{Soage González, José Carlos.}
#' \item{de Uña-Álvarez, Jacobo.}
#' \item{D. Hart, Jeffrey.}
#' }
#'
#' @references Cousido-Rocha, M., de Uña-Álvarez J., and Hart, J. (2018). A two-sample test for the equality of distributions for high-dimensional data. Preprint.
#'
#'
"_PACKAGE"
#> [1] "_PACKAGE"
