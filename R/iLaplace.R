#' iLaplace: A package for Approximating Multidimensional Integrals of Unimodal Functions
#'
#' This package gives (a parallel) implementation of the improved Laplace approximation for multivariate integrals of user-written unimodal functions (see Ruli et al. 2015). The method essentially approximates the target integral by the ratio of the (unnormalised) integrand and an approximation of its normalised version, both evaluated at the modal value. The normalised integrand is obtained through a sequential application of the Laplace approximation for marginal densities. Like the standard Laplace approximation, the improved Laplace approximation is a deterministic method which approximates intractable multidimensional integrals by (essentially) numerical optimisations. However, whith respect to the Laplace approximation, the improved Laplace involves scalar numerical integrations. Nevertheless, the improved Laplace approximation tends to be fast and extremely accurate, especially with skewed fat-tailed integrands.
#'
#' Currently the packages provides two functions \code{iLap} and \code{ilap2d} which perform approximation of \eqn{d}{d}-variate and bivariate integrals, respectively.
#'
#' @docType package
#' @name iLaplace-package
#' @author Erlis Ruli \email{erlisr@@yahoo.it}
#' @seealso \code{\link[iLaplace]{iLap}} and \code{\link[iLaplace]{iLap2d}} for more details and examples.
#' @references
#' Ruli E., Sartori N. & Ventura L. (2015)
#' Improved Laplace approximation for marignal likelihoods. \url{http://arxiv.org/abs/1502.06440}
NULL
