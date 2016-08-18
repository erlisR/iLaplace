##' @name iLap
##' @importFrom stats nlminb
##' @import doParallel
##' @import foreach
##' @import iterators
##' @import Rcpp
##' @import fastGHQuad
##' @useDynLib iLaplace
##'
##' @title Improved Laplace approximation for integrals of unimodal functions
##'
##' @description This function implements the improved Laplace approximation of Ruli et al. (2015) for multivariate integrals of user-written unimodal functions. See "Details" below for more information.
##' @usage iLap(fullOpt, ff, ff.gr, ff.hess, quad.data, ...)
##' @param fullOpt A list containing the minium (to be accesed via \code{fullOpt$par}), the value of the function at the minimum (to be accessed via \code{fullOpt$objective}) and the Hessian matrix at the minimum (to be accessed via \code{fullOpt$hessian}
##' @param ff The minus logarithm of the integrand function (the \code{h} function; see "Details").
##' @param ff.gr The gradient of \code{ff}, having the exact same arguments as  \code{ff}.
##' @param ff.hess The Hessian matrix of\code{ff}, having the exact same arguments as  \code{ff}.
##' @param quad.data Data for the Gaussian-Herimte quadratures; see "Details"
##' @param ... Additional arguments to be passed to \code{ff}, \code{ff.gr} and \code{ff.hess}
##' @return A double, the logarithm of the integral
##'
##' @details \code{iLap} approximates integrals of the type \deqn{I = \int_{x\in\mathcal{R}^d}\exp\{-h(x)\}\,dx}{I = \int\exp(-h(x)) dx} where \eqn{-h(\cdot)}{-h()} is a concave and unimodal function, with \eqn{x}{x} being \eqn{d}{d} dimensional real vector (\eqn{d>1}{d>1}). The approximation of \eqn{I} is obtained as the ratio between the unormalised kernel \eqn{-h(x)}{-h(x)} and an approximate density function \eqn{f(x)}{f(x)}, both evaluated at the modal value \eqn{x = \hat{x}}{x = \hat{x}}. The approximate density function \eqn{f(x)}{f(x)} is obtained by resorting to the Laplace approximation for marginal densities. The minimisations are performed with \code{\link[stats]{nlminb}} by suppling the gradient \code{ff.gr} and Hessian matrix {ff.hess} of \eqn{f(x)}{f(x)}. The normalisation of the univariate components is perforemd via Gaussian-Hermite quadratures as implemented in the function \code{\link[fastGHQuad]{aghQuad}}. The Gaussian-Quadrature data, to be provided via the argument \code{quad.data}, can be computed with the function \code{\link[fastGHQuad]{gaussHermiteData}} for a desired number of quadrature points. See "Examples" below.
##'
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}
##'
##' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite
##' Quadrature. \emph{Biometrika} \bold{81}, 624-629.
##' @examples
##'
##' # The negative integrand function in log
##' # is the negative log-density of the multivariate
##' # Student-t density centred at 0 with unit scale matrix
##' ff <- function(x, df) {
##'        d <- length(x)
##'        S <- diag(1, d, d)
##'        S.inv <- solve(S)
##'        Q <- colSums((S.inv %*% x) * x)
##'        logDet <- determinant(S)$modulus
##'        logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) +
##'        logDet) - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
##'        return(-logPDF)
##'        }
##'
##' # the gradient of ff
##' ff.gr <- function(x, df){
##'             m <- length(x)
##'             kr = 1 + crossprod(x,x)/df
##'             return((m+df)*x/(df*kr))
##'             }
##'
##' # the Hessian matrix of ff
##' ff.hess <- function(x, df) {
##' m <- length(x)
##' kr <- as.double(1 + crossprod(x,x)/df)
##' ll <- -(df+m)*2*tcrossprod(x,x)/(df*kr)^2.0
##' dd = (df+m)*(kr - 2*x^2/df)/(df*kr^2.0)
##' diag(ll) = dd;
##' return(ll)
##' }
##'
##' df = 5
##' dims <- 5:15

# improved Laplace and standard Laplace computation
##' normConts <- sapply(dims, function(mydim) {
##' opt <- nlminb(rep(1,mydim), ff, gradient = ff.gr, hessian = ff.hess, df = df)
##' opt$hessian <- ff.hess(opt$par, df = df);
##' quad.data = gaussHermiteData(50)
##' iLap <- iLap(opt, ff, ff.gr, ff.hess, quad.data = quad.data, df = df);
##' Lap <- mydim*log(2*pi)/2 - opt$objective - 0.5*determinant(opt$hessian)$mod;
##' return(c(iLap = iLap, Lap = Lap))
##' })

##' # plot the results
##' \dontrun{
##' plot(dims, normConts[1,], pch="*", cex = 1.6,
##'  ylim = c(-5, 0)) #improved Laplace
##' lines(dims, normConts[2,], type = "p", pch = "+") #standard Laplace
##' abline(h = 0) # the true value
##' }
##'
##'\dontrun{
##' ## See also the examples provided in the pacakge iLaplaceExamples, which is
##' ## an auxiliary R pacakge for iLaplace. To download it (be sure you have
##' ## the devtools package) run from R
##' ## devtools::install_github(erlisR/iLaplaceExamples)
##' ## or download the source at \url{https://github.com/erlisR/iLaplaceExamples}.
##'
##' }
##'
##'
##' @export
iLap <- function(fullOpt, ff, ff.gr, ff.hess, quad.data, ...)
{
  i = NULL

  m = length(fullOpt$par)
  obj = aux_quant(fullOpt$hessian, m)
  se = obj[[1]]
  # se = SEv(fullOpt$hessian, m)

  # fullOpt$ldblock =  ldetHessBlocks(fullOpt$hessian, m)
  fullOpt$ldblock =  obj[[m+1]]
  tmp = sapply(1:m, function(i) log(aghQuad(g = ila.densv,
                                            muHat = fullOpt$par[i],
                                            sigmaHat = se[i],
                                            rule =  quad.data,
                                            fullOpt = fullOpt,
                                            ff = ff,
                                            ff.gr = ff.gr,
                                            ff.hess = ff.hess,
                                            index = i,
                                            m = m, ...)))
  norm.const = sum(tmp)

  out = -m*0.5*log(2*pi) + 0.5*fullOpt$ldblock[1]

  return(-fullOpt$objective - out + norm.const)
}
