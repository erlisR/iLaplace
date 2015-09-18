# Given "q"-dimensional vector "par", "q"-vector "se", and scalars "lengthOut",
# "q" and "delta", this function creates for each of the "q" elements of
# "par" a grid of "lengthOut" values from "par-delta*se", to "par+delta*se". It
# returns a matrix of "lengthOut" times "q" values.
seqMat <- function(par, se, lengthOut, q, delta) {
  .Call('iLaplace_seqMat', PACKAGE = 'iLaplace', par, se, lengthOut, q, delta)
}

# Given the "dimMat" times "dimMat" matrix "hessMat" (being positve definite),
# this function computes the log-determinant of all the diagonal blocks,
# starting from the whole matrix. The function returns a "dimMat"-vector, where
# the first element is the log-determinant of hessMat, the second element is
# the log-determinant of hessMat[-1,-1] and so on, except the last element which
# is simply log(hessMat[dimMat,dimMat])
ldetHessBlocks <- function(hessMat, dimMat) {
  .Call('iLaplace_ldetHessBlocks', PACKAGE = 'iLaplace', hessMat, dimMat)
}

# Given the "dimMat" times "dimMat" matrix "hessMat" (being positve definite),
# this function computes the square root of the diagonal of inverse of blocks
# of "hessMat". The function returns a "dimMat" vector with the first element
# being sqrt(diag(solve(hessMat)))[1], the second element is sqrt(diag(solve(hessMat[-1,-1])))[1],
# and so on, except fo the last element which is sqrt(1/hessMat[dimMat,dimMat]).
SEv <- function(hessMat, dimMat) {
  .Call('iLaplace_SEv', PACKAGE = 'iLaplace', hessMat, dimMat)
}
