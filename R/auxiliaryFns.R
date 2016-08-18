# Given "q"-dimensional vector "par", "q"-vector "se", and scalars "lengthOut",
# "q" and "delta", this function creates for each of the "q" elements of
# "par" a grid of "lengthOut" values from "par-delta*se", to "par+delta*se". It
# returns a matrix of "lengthOut" times "q" values.
# seqMat <- function(par, se, lengthOut, q, delta) {
#   .Call('iLaplace_seqMat', PACKAGE = 'iLaplace', par, se, lengthOut, q, delta)
# }

# Given the "dimMat" times "dimMat" matrix "hessMat" (being positve definite),
# this function computes the log-determinant of all the diagonal blocks,
# starting from the whole matrix. The function returns a "dimMat"-vector, where
# the first element is the log-determinant of hessMat, the second element is
# the log-determinant of hessMat[-1,-1] and so on, except the last element which
# is simply log(hessMat[dimMat,dimMat])
# ldetHessBlocks <- function(hessMat, dimMat) {
#   .Call('iLaplace_ldetHessBlocks', PACKAGE = 'iLaplace', hessMat, dimMat)
# }

# Given the "dimMat" times "dimMat" matrix "hessMat" (being positve definite),
# this function computes the square root of the diagonal of inverse of blocks
# of "hessMat". The function returns a "dimMat" vector with the first element
# being sqrt(diag(solve(hessMat)))[1], the second element is sqrt(diag(solve(hessMat[-1,-1])))[1],
# and so on, except fo the last element which is sqrt(1/hessMat[dimMat,dimMat]).
# SEv <- function(hessMat, dimMat) {
#   .Call('iLaplace_SEv', PACKAGE = 'iLaplace', hessMat, dimMat)
# }
#
# fnblocks <- function(hessMat, dimMat) {
#   .Call('iLaplace_fnblocks', PACKAGE = 'iLaplace', hessMat, dimMat)
# }

# objfun <- function(argv, fullOpt, ff, ff.gr, ff.hess, m, lo, up, control, ...){
#   xv <- argv[1:control$sp.points]
#   index <- argv[control$sp.points + 1]
#
#   fun.val <- sapply(xv, ila.dens, fullOpt = fullOpt, ff = ff,
#                     ff.gr = ff.gr, ff.hess = ff.hess, index = index,
#                     m = m, ...)
#
#   den.sp <- splinefun(x = xv, y = fun.val)
#   #   plot(cond.sp, lo[index], up[index])
#
#   nc.den  <-  integrate(den.sp, lower = lo[index], upper = up[index])$value
#
#   log.den  <-  log(ila.dens(fullOpt$par[index], fullOpt = fullOpt, ff = ff,
#                             ff.gr = ff.gr, ff.hess = ff.hess, index = index,
#                             m = m, ...))
#
#   return(log.den - log(nc.den))
# }

ila.dens <- function(pp, fullOpt, ff, ff.gr, ff.hess, index, m, ...) {

  if (index == 1) {
    ans <- marg(x = pp, fullOpt = fullOpt, ff = ff, ff.gr = ff.gr, ff.hess = ff.hess, ...)
  } else {
    if (index < m) {
      ans <- middleConds(x = pp, fullOpt = fullOpt, ff = ff, ff.gr = ff.gr,
                         ff.hess = ff.hess, index = index, m = m, ...)
    } else {
      ans <- lastCond(x = pp, fullOpt = fullOpt, ff = ff, m = m, ...)
    }
  }
  return(ans)
}

ila.densv <- Vectorize(ila.dens, "pp")

ila.densCW <- function(pp, fullOpt, ff, ff.gr, ff.hess, index, m, ...) {

  if (index == 1) {
    ans <- margCW(x = pp, fullOpt = fullOpt, ff = ff, ff.gr = ff.gr, ff.hess = ff.hess, ...)
  } else {
    if (index < m) {
      ans <- middleCondsCW(x = pp, fullOpt = fullOpt, ff = ff, ff.gr = ff.gr,
                         ff.hess = ff.hess, index = index, m = m, ...)
    } else {
      ans <- lastCond(x = pp, fullOpt = fullOpt, ff = ff, m = m, ...)
    }
  }
  return(ans)
}

ila.densCWv <- Vectorize(ila.densCW, "pp")

marg <- function(x, fullOpt, ff, ff.gr, ff.hess, ...) {

  out = 0.0

  try({
    tmp = function(xx, ...) ff(c(x, xx), ...)
    gr.tmp = function(xx, ...) ff.gr(c(x, xx), ...)[-1]
    hess.tmp = function(xx, ...) ff.hess(c(x, xx), ...)[-1, -1]
    optTmp = nlminb(start = fullOpt$par[-1], objective = tmp,
                    gradient = gr.tmp, hessian = hess.tmp, ...)
    optTmp$hessian = ff.hess(c(x, optTmp$par), ...)[-1, -1]
    ldetc = determinant(optTmp$hessian)$mod
    out = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[1] -
                                        ldetc) - optTmp$obj + fullOpt$obj
    out <- exp(out)

    if(!is.finite(out))
      out = 0.0
  })

  return(out)
}

margCW <- function(x, fullOpt, ff, ff.gr, ff.hess, ...) {

  out = 0.0
  try({
    cstOpt = list()
    # tmp = function(xx, ...) ff(c(x, xx), ...)
    # gr.tmp = function(xx, ...) ff.gr(c(x, xx), ...)[-1]
    # hess.tmp = function(xx, ...) ff.hess(c(x, xx), ...)[-1, -1]
    # optTmp = nlminb(start = fullOpt$par[-1], objective = tmp,
    #                 gradient = gr.tmp, hessian = hess.tmp, ...)
    cstOpt$par = fullOpt$par[-1] + fullOpt$delta[[1]]*(fullOpt$par[1]-x);

    cstOpt$objective = ff(c(x, cstOpt$par), ...)
    cstOpt$hessian = ff.hess(c(x, cstOpt$par), ...)[-1, -1]
    ldetc = determinant(cstOpt$hessian)$mod

    out = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[1] -
                                        ldetc) - cstOpt$obj + fullOpt$obj
    out <- exp(out)

    if(!is.finite(out))
      out = 0.0
  })

  return(out)
}

middleConds <- function(x, fullOpt, ff, ff.gr, ff.hess, index, m, ...) {

  out = 0.0

  try({
    no = 1:index
    tmp = function(xx, ...) ff(c(fullOpt$par[1:(index - 1)], x,
                                 xx), ...)
    gr.tmp = function(xx, ...) ff.gr(c(fullOpt$par[1:(index -
                                                        1)], x, xx), ...)[-no]
    hess.tmp = function(xx, ...) ff.hess(c(fullOpt$par[1:(index - 1)], x, xx), ...)[-no, -no]
    if (index == (m - 1)) {
      tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp,
                      gradient = gr.tmp, hessian = NULL, ...)
      tmpOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
                                 x, tmpOpt$par), ...)[-no, -no]
      ldetc = log(tmpOpt$hessian)
    }
    else {
      tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp,
                      gradient = gr.tmp, hessian = hess.tmp, ...)
      tmpOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
                                 x, tmpOpt$par), ...)[-no, -no]

      ldetc = as.double(determinant(tmpOpt$hessian)$mod)
    }
    ans = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[index] -
                                        ldetc) - tmpOpt$obj + fullOpt$obj
    out = exp(ans)

    if(!is.finite(out)){
      out = 0.0
    }

  })

  return(out)
}

middleCondsCW <- function(x, fullOpt, ff, ff.gr, ff.hess, index, m, ...) {
# handles conditionals for 2 <= index <= m-1
  out = 0.0

  try({
    cstOpt = list()
    no = 1:index
    # tmp = function(xx, ...) ff(c(fullOpt$par[1:(index - 1)], x,
    #                              xx), ...)
    # gr.tmp = function(xx, ...) ff.gr(c(fullOpt$par[1:(index -
    #                                                     1)], x, xx), ...)[-no]

    # hess.tmp = function(xx, ...) ff.hess(c(fullOpt$par[1:(index - 1)], x, xx), ...)[-no, -no]

    if (index == (m - 1)) {
      # tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp,
      #                 gradient = gr.tmp, hessian = NULL, ...)
      # tmpOpt = fullOpt$par[-no] + fullOpt$hessian[m,index]*(fullOpt$par[index]-x)/fullOpt$hessian[m,m]
      cstOpt$par = fullOpt$par[-no] + fullOpt$delta[[index]]*(fullOpt$par[index]-x);
      cstOpt$objective = ff(c(fullOpt$par[1:(index - 1)],
                              x, cstOpt$par), ...)
      cstOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
                                 x, cstOpt$par), ...)[-no, -no]
      ldetc = log(cstOpt$hessian)

    } else {
      # optTmp$par = lamb_psi(x, fullOpt$par[-no], fullOpt$par[index],
      #                       fullOpt$invHess[-no,-no], fullOpt$hessian[-no,no],
      #                       length(fullOpt$par[-no]))
      cstOpt$par = fullOpt$par[-no] + fullOpt$delta[[index]]*(fullOpt$par[index]-x);

      cstOpt$objective = ff(c(fullOpt$par[1:(index - 1)],
                              x, cstOpt$par), ...)
      cstOpt$hessian = ff.hess(c(fullOpt$par[1:(index - 1)],
                                 x, cstOpt$par), ...)[-no, -no]

      ldetc = as.double(determinant(cstOpt$hessian)$mod)
    }

    ans = -0.5 * log(2 * pi) + 0.5 * (fullOpt$ldblock[index] -
                                        ldetc) - cstOpt$obj + fullOpt$obj
    out = exp(ans[1])

    if(!is.finite(out)){
      out = 0.0
    }

  })

  return(out)
}

lastCond <- function(x, fullOpt, ff, m, ...) {

  out <- 0.0

  try({
    tmp = function(x) ff(c(fullOpt$par[c(1:(m - 1))], x), ...)
    out = -0.5 * log(2 * pi) + 0.5 * fullOpt$ldblock[m] - tmp(x) + fullOpt$obj

    out <- exp(out)

    if (!is.finite(out))
      out <- 0.0
    })

  return(out)
}
