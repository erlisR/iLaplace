##' @name iLap2d
##' @title Improved Laplace approximation for bivariate integrals of unimodal functions
##' @description This function is similar to \code{iLap} except that it handles only bivariate integrals of user-written unimodal functions.
##' @usage iLap2d(fullOpt, ff, ff.gr, ff.hess, quad.data, ...)
##'
##' @param fullOpt A list containing the minium (to be accesed via \code{fullOpt$par}), the value of the function at the minimum (to be accessed via \code{fullOpt$objective}) and the Hessian matrix at the minimum (to be accessed via \code{fullOpt$hessian}
##' @param ff The minus logarithm of the integrand function (the \code{h} function, see \code{iLap} for further details).
##' @param ff.gr The gradient of \code{ff}, having the exact same arguments as  \code{ff}
##' @param ff.hess The Hessian matrix of\code{ff}, having the exact same arguments as  \code{ff}
##' @param quad.data Data for the Gaussian-Herimte quadratures; see "Details"
##' @param ... Additional arguments to be passed to \code{ff}, \code{ff.gr} and \code{ff.hess}
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}
##'
##' @examples
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
##' dgf = 5
##' opt <- nlminb(rep(1,2), ff, gradient = ff.gr, hessian = ff.hess, df = dgf)
##' opt$hessian <- ff.hess(opt$par, df = dgf);
##' quad.data = gaussHermiteData(50)
##'
##' # The improved Laplace approximation (the truth equals 0.0)
##' iLap <- iLap2d(fullOpt = opt, ff = ff, ff.gr = ff.gr,
##'                ff.hess = ff.hess, quad.data = quad.data,
##'                df = dgf)
##' # The standard Laplace approximation (the truth equals 0.0)
##' Lap <- log(2*pi) - opt$objective - 0.5*determinant(opt$hessian)$mod;
##'
##' @return a double, the logarithm of the integral
##' @export
iLap2d <- function(fullOpt, ff, ff.gr, ff.hess, quad.data, ...)
{
  m = 2
  obj = aux_quant(fullOpt$hessian, m)
  se = obj[[1]]
  fullOpt$ldblock =  obj[[m+1]]

  # the marginal
  marg2 = function(x1, ...) {
    out = 0.0;
    try({
      tmp = function(x) ff(c(x1, x), ...)
      gr.tmp = function(x) ff.gr(c(x1, x), ...)[-1]
      hes.tmp = function(x) ff.hess(c(x1,x), ...)[-1,-1]
      tmpOpt = nlminb(fullOpt$par[2], tmp, gradient = gr.tmp)
      tmpOpt$hessian = hes.tmp(tmpOpt$par)
      out = exp(-0.5*log(2*pi) + 0.5*(fullOpt$ldblock[1] -
                                        log(abs(tmpOpt$hessian))) -
                  tmpOpt$obj + fullOpt$obj)

      if(!is.finite(out))
        out = 0.0
    })
    return(out)
  }

  marg2v = Vectorize(marg2, vectorize.args = "x1")

  nc.marg  = log(aghQuad(g = marg2v,
                         muHat = fullOpt$par[1],
                         sigmaHat = se[1],
                         rule =  quad.data,
                         ...))


  # the conditional
  cond2 <- function(x, ...) {

    out <- 0.0

    try({
      tmp = function(x) ff(c(fullOpt$par[c(1:(m - 1))], x), ...)
      out = -0.5 * log(2 * pi) + 0.5 * fullOpt$ldblock[m] - tmp(x) + fullOpt$obj

      out = exp(out)

      if (!is.finite(out))
        out <- 0.0
    })

    return(out)
  }
  cond2v <- Vectorize(cond2, vectorize.args = "x")

  nc.cond  = log(aghQuad(g = cond2v,
                     muHat = fullOpt$par[2],
                     sigmaHat = se[2],
                     rule =  quad.data,
                     ...))

  out = -m*0.5*log(2*pi) + 0.5*fullOpt$ldblock[1]
  norm.const = nc.marg + nc.cond

  return(-fullOpt$objective - out + norm.const)
}
