##' @name iLap_par

##' @title Improved Laplace approximation in parallel
##'
##' @description Does the same as \code{iLap} but in parallel
##' @usage iLap_par(fullOpt, ff, ff.gr, ff.hess, quad.data,
##'      control = list(n.cores = 1), ...)
##' @param fullOpt A list containing the minium (to be accesed via \code{fullOpt$par}), the value of the function at the minimum (to be accessed via \code{fullOpt$objective}) and the Hessian matrix at the minimum (to be accessed via \code{fullOpt$hessian}
##' @param ff The minus logarithm of the integrand function (the \code{h} function; see "Details").
##' @param ff.gr The gradient of \code{ff}, having the exact same arguments as  \code{ff}.
##' @param ff.hess The Hessian matrix of\code{ff}, having the exact same arguments as  \code{ff}.
##' @param quad.data Data for the Gaussian-Herimte quadratures; see "Details"
##' @param control A named list of control parameters with elements \code{n.cores} which sets the number of cores to be used for the parallel computiations. See "Details" for more information.
##' @param ... Additional arguments to be passed to \code{ff}, \code{ff.gr} and \code{ff.hess}
##' @return A double, the logarithm of the integral
##'
##' @details See \code{iLap}.
##'
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}
##'
##' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite
##' Quadrature. \emph{Biometrika} \bold{81}, 624-629.
##'
##' Cox, D.R and Wermuth, W. (1990). An approximation to maximum
##' likelihood estimates in reduced models. \emph{Biometrika} \bold{77}, 747-761
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
##' iLap <- iLap_par(opt, ff, ff.gr, ff.hess, quad.data = quad.data, df = df);
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
iLap_par <- function(fullOpt, ff, ff.gr, ff.hess, quad.data,
                 control = list(n.cores = 1),
                 ...)
{
  i = NULL

  m = length(fullOpt$par)
  obj = aux_quant(fullOpt$hessian, m)
  se = obj[[1]]
  # se = SEv(fullOpt$hessian, m)

  # fullOpt$ldblock =  ldetHessBlocks(fullOpt$hessian, m)
  fullOpt$ldblock =  obj[[m+1]]
  registerDoParallel(cores = control$n.cores)

  norm.const <- foreach(i = 1:m, .combine = sum) %dopar% {
    log(aghQuad(g = ila.densv,
                muHat = fullOpt$par[i],
                sigmaHat = se[i],
                rule =  quad.data,
                fullOpt = fullOpt,
                ff = ff,
                ff.gr = ff.gr,
                ff.hess = ff.hess,
                index = i,
                m = m, ...))
    }

  out = -m*0.5*log(2*pi) + 0.5*fullOpt$ldblock[1]

  return(-fullOpt$objective - out + norm.const)
}
