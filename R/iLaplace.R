##' @name iLaplace
##' @import parallel
##' @useDynLib iLaplace
##'
##' @title Improved Laplace approximation (using analytical gradient and Hessian)
##'
##' @description This function computes the improved Laplace approximation of Ruli et al. (2015) for multivariate integrals of user-written functions. See Details below for more information. For practical examples on the use of this package please refer to the \code{iLaplaceExamples} package on \url{https://github.com/erlisR/iLaplaceExamples}.
##' @usage iLaplace(fullOpt, ff, ff.gr, ff.hess,
##'          control = list(sp.points = 100, delta = 13, n.cores = detectCores()-1),
##'          clEvalQ = NULL, ...)
##' @param fullOpt A list containing the minium (to be accesed via \code{fullOpt$par}), the value of the function at the minimum (to be accessed via \code{fullOpt$objective}) and the Hessian matrix at the minimum (to be accessed via \code{fullOpt$hessian}
##' @param ff The minus logarithm of the integrand function (the \code{h} function, see Details).
##' @param ff.gr The gradient of \code{ff}, having the exact same arguments as  \code{ff}
##' @param ff.hess The Hessian matrix of\code{ff}, having the exact same arguments as  \code{ff}
##' @param control A named list of control parameters with elements \code{sp.points}, \code{delta} and \code{n.cores}. \code{sp.points} sets the number points for the spline evaluations; \code{delta} controls the length of the inteval of integration; \code{n.cores} sets the number of cores to be used for the parallel computiations. This is available for Linux/Unix systems only as it uses forking. See Details for more information.
##' @param clEvalQ Other useful R packages to be passed to the (internally crated) clusters for the parallel computations. If required, provide package names as a vector of strings.
##' @param ... Additional arguments to be passed to \code{ff}, \code{ff.gr} and \code{ff.hess}
##' @return double, the logarithm of the integral
##' @details \code{iLaplace} approximates integrals of the form \deqn{I = \int_{x\in\mathcal{R}^d}\exp\{-h(x)\}\,dx}{I = \int\exp(-h(x)) dx} where \eqn{-h(\cdot)}{-h()} is a concave and unimodal function, with \eqn{x}{x} being \eqn{d}{d} dimensional real vector (\eqn{d>1}{d>1}). The approximation of \eqn{I} is obtained as the ratio between the unormalised kernel \eqn{-h(x)}{-h(x)} and an approximate density function \eqn{f(x)}{f(x)}, both evaluated at the modal value \eqn{x = \hat{x}}{x = \hat{x}}. The approximate density function \eqn{f(x)}{f(x)} is obtained by resorting to the Laplace approximation for marginal densities. The normalisation of the univariate components is done as follows. First, for each of the \eqn{d}{d} dimensions, a suitable grid is fixed in the neighbourhood of \eqn{\hat x}{\hat x}. To be sure all the region with non-negligible mass is considered, the extreemes of the neighbourhood are set to \eqn{(\hat x - \delta se,\hat x + \delta se)}{(\hat x - delta*se, \hat x - delta*se)} and the grid is made of \code{sp.points} equally spaced values. Here, \eqn{se}{se} is a suitable "profile" standard error that is computed from the Hessian matrix evaluated at the modal value. \code{delta} is fixed to 13 by default. Finally, the scalar Laplace approximations for marginal densities are evaluated over the grids and the densities are reconstrcuted via the \code{\link[stats]{splinefun}} function. Lastly the splines are then normalized with the function \code{\link[stats]{integrate}} over the region constructed as above.
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2015)
##' Improved Laplace approximation for marignal likelihoods.
##' \url{http://arxiv.org/abs/1502.06440}
##' @examples
##'\dontrun{
##'
##' ## See the examples provided in the pacakge iLaplaceExamples, which is
##' ## an auxiliary R pacakge for iLaplace. To download it run
##' ## devtools::install_github(erlisR/iLaplaceExamples) in R or go
##' ## to \url{https://github.com/erlisR/iLaplaceExamples}.
##'
##' }
##'
##'
##' @export
# iLaplace approximation for d-variate integral, with user-supplied objective function
iLaplace <- function(fullOpt, ff, ff.gr, ff.hess,
                     control = list(sp.points = 100, delta = 13, n.cores = detectCores() - 1),
                     clEvalQ = NULL, ...)
{
  m = length(fullOpt$par)
  se = SEv(fullOpt$hessian, m)
  oo = seqMat(par = fullOpt$par, se, lengthOut = control$sp.points, q = m, delta = control$delta)
  par.val <- matrix(NA, control$sp.points + 1, m)
  par.val[1:control$sp.points,] <- oo$parVal
  par.val[control$sp.points + 1,] <- 1:m
  fullOpt$ldblock =  ldetHessBlocks(fullOpt$hessian, m)

  objfun <- function(argv, ...){
    lo = oo$lo
    up = oo$up
    xv <- argv[1:control$sp.points]
    index <- argv[control$sp.points + 1]
    ila.dens <- function(pp, index) {
      marg = function(x) {
        out <- tryCatch({
          tmp = function(xx) ff(c(x, xx), ...)
          gr.tmp = function(xx) ff.gr(c(x, xx), ...)[-1]
          hess.tmp = function(xx) ff.hess(c(x, xx), ...)[-1, -1]
          optTmp = nlminb(start = fullOpt$par[-1], objective = tmp,
                          gradient = gr.tmp, hessian = hess.tmp)
          optTmp$hessian = ff.hess(c(x, optTmp$par), ...)[-1, -1]
          ldetc = determinant(optTmp$hessian)$mod
          out = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[1] -
                                              ldetc) - optTmp$obj + fullOpt$obj
          exp(as.double(out))
        },
        error=function(cond) {
          return(NA)
        },
        warning=function(cond) {
          return(NULL)
        },
        finally={

        }
        )
        return(out)
      }

      middleConds = function(x, index) {
        out=tryCatch({
          no = c(1:index)
          tmp = function(xx) ff(c(fullOpt$par[1:(index - 1)], x,
                                  xx), ...)
          gr.tmp = function(xx) ff.gr(c(fullOpt$par[1:(index -
                                                         1)], x, xx), ...)[-no]
          hess.tmp = function(xx) ff.hess(c(fullOpt$par[1:(index - 1)], x, xx), ...)[-no, -no]
          if (index == (m - 1)) {
            tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp,
                            gradient = gr.tmp)
            tmpOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
                                       x, tmpOpt$par), ...)[-no, -no]
            ldetc = log(tmpOpt$hessian)
          }
          else {
            tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp,
                            gradient = gr.tmp, hessian = hess.tmp)
            tmpOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
                                       x, tmpOpt$par), ...)[-no, -no]

            ldetc = as.double(determinant(tmpOpt$hessian)$mod)
          }
          ans = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[index] -
                                              ldetc) - tmpOpt$obj + fullOpt$obj
          exp(ans)
        },
        error=function(cond) {
          return(0.0)
        },
        warning=function(cond) {
          return(NULL)
        },
        finally={
        }
        )
        return(out)
      }

      lastCond = function(x) {
        out=tryCatch({
          tmp = function(x) ff(c(fullOpt$par[c(1:(m - 1))], x), ...)
          out = -0.5 * log(2 * pi) + 0.5 * log(fullOpt$hessian[m, m]) - tmp(x) + fullOpt$obj
          exp(out)
        },
        error=function(cond) {
          return(0.0)
        },
        warning=function(cond) {
          return(NULL)
        },
        finally={
        }
        )
        return(out)
      }

      if(index==1) {
        out = marg(pp)
      } else {
        if(index>1 && index < m){
          out = middleConds(pp, index)
        } else {
          out = lastCond(pp)
        }
      }
      return(out)
    }
    cfun.val <- sapply(xv, function(x) ila.dens(x, index))
    cond.sp <- splinefun(x = xv, y = cfun.val)
    #   plot(cond.sp, lo[index], up[index])
    nc.cond  <-  integrate(cond.sp, lower = lo[index], upper = up[index])$value
    log.den  <-  log(ila.dens(fullOpt$par[index], index))
    return(log.den - log(nc.cond))
  }

  if (control$n.cores < 2) {
    message("#----------------------------------------\n cores < 2 -> computing the integral serially...\n#-----------------------------------------")
    ilaf <- apply(X = par.val, MARGIN = 2, objfun, ... )

  } else {

    cl <- makeCluster(control$n.cores, type = "SOCK")

    ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }

    clusterExport(cl = cl, ls(), envir = environment())

    clusterEvalQ(cl = cl, {
      ipak <- function(pkg){
      new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
      if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
      sapply(pkg, require, character.only = TRUE)
    };
      ipak(clEvalQ);})

    ilaf <- parApply(cl = cl, X = par.val, MARGIN = 2, objfun, ... )

    stopCluster(cl = cl)
  }

  return(-fullOpt$obj - sum(ilaf))
}

# iLap_an2 <- function(fullOpt, ff, ff.gr, ff.hess,
#                      control = list(sp.points=100,
#                                     delta=13,
#                                     n.cores = detectCores() - 1))
# {
#   m = length(fullOpt$par)
#   se = SEv(fullOpt$hessian, m)
#   oo = seqMat(par=fullOpt$par, se, lengthOut=control$sp.points, q=m, delta=control$delta)
#   lo = oo$lo
#   up = oo$up
#   par.val <- oo$parVal
#   fullOpt$ldblock =  ldetHessBlocks(fullOpt$hessian, m)
#
#
#   ila.dens <- function(pp, index) {
#     marg = function(x) {
#       out <- tryCatch({
#         tmp = function(xx) ff(c(x, xx))
#         gr.tmp = function(xx) ff.gr(c(x, xx))[-1]
#         hess.tmp = function(xx) ff.hess(c(x, xx))[-1, -1]
#         optTmp = nlminb(start = fullOpt$par[-1], objective = tmp,
#                         gradient = gr.tmp)
#         optTmp$hessian = ff.hess(c(x, optTmp$par))[-1, -1]
#         ldetc = determinant(optTmp$hessian)$mod
#         out = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[1] -
#                                             ldetc) - optTmp$obj + fullOpt$obj
#         exp(as.double(out))
#       },
#       error=function(cond) {
#         return(NA)
#       },
#       warning=function(cond) {
#         return(NULL)
#       },
#       finally={
#
#       }
#       )
#       return(out)
#     }
#
#     middleConds = function(x, index) {
#       out=tryCatch({
#         no = c(1:index)
#         tmp = function(xx) ff(c(fullOpt$par[1:(index - 1)], x,
#                                 xx))
#         gr.tmp = function(xx) ff.gr(c(fullOpt$par[1:(index -
#                                                        1)], x, xx))[-no]
#         hess.tmp = function(xx) ff.hess(c(fullOpt$par[1:(index - 1)], x, xx))[-no, -no]
#         if (index == (m - 1)) {
#           tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp,
#                           gradient = gr.tmp)
#           tmpOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
#                                      x, tmpOpt$par))[-no, -no]
#           ldetc = log(tmpOpt$hessian)
#         }
#         else {
#           tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp,
#                           gradient = gr.tmp)
#           tmpOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
#                                      x, tmpOpt$par))[-no, -no]
#
#           ldetc = as.double(determinant(tmpOpt$hessian)$mod)
#         }
#         ans = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[index] -
#                                             ldetc) - tmpOpt$obj + fullOpt$obj
#         exp(ans)
#       },
#       error=function(cond) {
#         return(0.0)
#       },
#       warning=function(cond) {
#         return(NULL)
#       },
#       finally={
#       }
#       )
#       return(out)
#     }
#
#     lastCond = function(x) {
#       out=tryCatch({
#         tmp = function(x) ff(c(fullOpt$par[c(1:(m - 1))], x))
#         out = -0.5 * log(2 * pi) + 0.5 * log(fullOpt$hessian[m, m]) - tmp(x) + fullOpt$obj
#         exp(out)
#       },
#       error=function(cond) {
#         return(0.0)
#       },
#       warning=function(cond) {
#         return(NULL)
#       },
#       finally={
#       }
#       )
#       return(out)
#     }
#
#     if(index==1) {
#       out = marg(pp)
#     } else {
#       if(index>1 && index < m){
#         out = middleConds(pp, index)
#       } else {
#         out = lastCond(pp)
#       }
#     }
#     return(out)
#   }
#
#   ila.hat <- function(xv, index){
#     cfun.val <- sapply(xv, function(x) ila.dens(x, index))
#     cond.sp <- splinefun(x = xv, y = cfun.val)
#     #   plot(cond.sp, lo[index], up[index])
#     nc.cond  <-  integrate(cond.sp, lower = lo[index], upper = up[index])$value
#     log.den  <-  log(ila.dens(fullOpt$par[index], index))
#     return(log.den-log(nc.cond))
#   }
#
#   #   plot(Vectorize(function(x) ila.hat(x, 1), "x"), lo[1], up[1])
#
#   ilaf <- 0.0
#   cl <- makeCluster(control$n.cores, type = "SNOW")
#   registerDoSNOW(cl)
#   i = 1
#   ilaf = sum(ilaf, foreach(i = 1:m, .combine = sum) %dopar% {
#     ilaf = ila.hat(par.val[,i], index = i)
#   })
#   return(-fullOpt$obj - ilaf)
# }
