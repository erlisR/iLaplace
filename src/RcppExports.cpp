// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// nlpost_gomp
double nlpost_gomp(NumericVector param, NumericVector data);
RcppExport SEXP iLaplace_nlpost_gomp(SEXP paramSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    __result = Rcpp::wrap(nlpost_gomp(param, data));
    return __result;
END_RCPP
}
// grad_gomp
NumericVector grad_gomp(NumericVector param, NumericVector data);
RcppExport SEXP iLaplace_grad_gomp(SEXP paramSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    __result = Rcpp::wrap(grad_gomp(param, data));
    return __result;
END_RCPP
}
// hess_gomp
NumericMatrix hess_gomp(NumericVector param, NumericVector data);
RcppExport SEXP iLaplace_hess_gomp(SEXP paramSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    __result = Rcpp::wrap(hess_gomp(param, data));
    return __result;
END_RCPP
}
// nlDen_mvt
double nlDen_mvt(arma::vec x, int m, double df);
RcppExport SEXP iLaplace_nlDen_mvt(SEXP xSEXP, SEXP mSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(nlDen_mvt(x, m, df));
    return __result;
END_RCPP
}
// nlDen_mvskt
double nlDen_mvskt(arma::vec x, double a, double c, int m, double df);
RcppExport SEXP iLaplace_nlDen_mvskt(SEXP xSEXP, SEXP aSEXP, SEXP cSEXP, SEXP mSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(nlDen_mvskt(x, a, c, m, df));
    return __result;
END_RCPP
}
// grad_mvt
arma::vec grad_mvt(arma::vec x, int m, double df);
RcppExport SEXP iLaplace_grad_mvt(SEXP xSEXP, SEXP mSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(grad_mvt(x, m, df));
    return __result;
END_RCPP
}
// grad_mvskt
arma::vec grad_mvskt(arma::vec x, double a, double c, int m, double df);
RcppExport SEXP iLaplace_grad_mvskt(SEXP xSEXP, SEXP aSEXP, SEXP cSEXP, SEXP mSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(grad_mvskt(x, a, c, m, df));
    return __result;
END_RCPP
}
// hess_mvt
arma::mat hess_mvt(arma::vec x, int m, double df);
RcppExport SEXP iLaplace_hess_mvt(SEXP xSEXP, SEXP mSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(hess_mvt(x, m, df));
    return __result;
END_RCPP
}
// hess_mvskt
arma::mat hess_mvskt(arma::vec x, double a, double c, int m, double df);
RcppExport SEXP iLaplace_hess_mvskt(SEXP xSEXP, SEXP aSEXP, SEXP cSEXP, SEXP mSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(hess_mvskt(x, a, c, m, df));
    return __result;
END_RCPP
}
// dmvt
double dmvt(arma::vec x, arma::vec mu, arma::mat S, int p, double df, bool lg);
RcppExport SEXP iLaplace_dmvt(SEXP xSEXP, SEXP muSEXP, SEXP SSEXP, SEXP pSEXP, SEXP dfSEXP, SEXP lgSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    __result = Rcpp::wrap(dmvt(x, mu, S, p, df, lg));
    return __result;
END_RCPP
}
// dhalfCauchy
double dhalfCauchy(double x, double scale, bool lg);
RcppExport SEXP iLaplace_dhalfCauchy(SEXP xSEXP, SEXP scaleSEXP, SEXP lgSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    __result = Rcpp::wrap(dhalfCauchy(x, scale, lg));
    return __result;
END_RCPP
}
// rmvnorm
arma::vec rmvnorm(arma::vec mu, arma::mat S, int p);
RcppExport SEXP iLaplace_rmvnorm(SEXP muSEXP, SEXP SSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    __result = Rcpp::wrap(rmvnorm(mu, S, p));
    return __result;
END_RCPP
}
// rmvt
arma::vec rmvt(arma::vec mu, arma::mat S, int p, double df);
RcppExport SEXP iLaplace_rmvt(SEXP muSEXP, SEXP SSEXP, SEXP pSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(rmvt(mu, S, p, df));
    return __result;
END_RCPP
}
// prJeff
double prJeff(double nu, bool lg);
RcppExport SEXP iLaplace_prJeff(SEXP nuSEXP, SEXP lgSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    __result = Rcpp::wrap(prJeff(nu, lg));
    return __result;
END_RCPP
}
// grldmvt
arma::vec grldmvt(arma::vec x, arma::vec mu, arma::mat S, int p, double df);
RcppExport SEXP iLaplace_grldmvt(SEXP xSEXP, SEXP muSEXP, SEXP SSEXP, SEXP pSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(grldmvt(x, mu, S, p, df));
    return __result;
END_RCPP
}
// hessldmvt
arma::mat hessldmvt(arma::vec x, arma::vec mu, arma::mat S, int p, double df);
RcppExport SEXP iLaplace_hessldmvt(SEXP xSEXP, SEXP muSEXP, SEXP SSEXP, SEXP pSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    __result = Rcpp::wrap(hessldmvt(x, mu, S, p, df));
    return __result;
END_RCPP
}
// nlpost_lub
double nlpost_lub(arma::vec beta, double lsig, arma::vec y, arma::vec x1, arma::vec x2, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_nlpost_lub(SEXP betaSEXP, SEXP lsigSEXP, SEXP ySEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(nlpost_lub(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// nlpostT_lub
double nlpostT_lub(arma::vec beta, double lsig, double lnu, arma::vec y, arma::vec x1, arma::vec x2, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_nlpostT_lub(SEXP betaSEXP, SEXP lsigSEXP, SEXP lnuSEXP, SEXP ySEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< double >::type lnu(lnuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(nlpostT_lub(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// grad_lub
arma::vec grad_lub(arma::vec beta, double lsig, arma::vec y, arma::vec x1, arma::vec x2, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_grad_lub(SEXP betaSEXP, SEXP lsigSEXP, SEXP ySEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(grad_lub(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// gradT_lub
arma::vec gradT_lub(arma::vec beta, double lsig, double lnu, arma::vec y, arma::vec x1, arma::vec x2, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_gradT_lub(SEXP betaSEXP, SEXP lsigSEXP, SEXP lnuSEXP, SEXP ySEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< double >::type lnu(lnuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(gradT_lub(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// hess_lub
arma::mat hess_lub(arma::vec beta, double lsig, arma::vec y, arma::vec x1, arma::vec x2, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_hess_lub(SEXP betaSEXP, SEXP lsigSEXP, SEXP ySEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(hess_lub(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// hessT_lub
arma::mat hessT_lub(arma::vec beta, double lsig, double lnu, arma::vec y, arma::vec x1, arma::vec x2, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_hessT_lub(SEXP betaSEXP, SEXP lsigSEXP, SEXP lnuSEXP, SEXP ySEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< double >::type lnu(lnuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(hessT_lub(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// MCMCmetrop_cpp
arma::mat MCMCmetrop_cpp(Function lfun, arma::vec start, arma::mat V, int mcmc, int burnin, double df, int verbose, Environment myframe);
RcppExport SEXP iLaplace_MCMCmetrop_cpp(SEXP lfunSEXP, SEXP startSEXP, SEXP VSEXP, SEXP mcmcSEXP, SEXP burninSEXP, SEXP dfSEXP, SEXP verboseSEXP, SEXP myframeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Function >::type lfun(lfunSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type mcmc(mcmcSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Environment >::type myframe(myframeSEXP);
    __result = Rcpp::wrap(MCMCmetrop_cpp(lfun, start, V, mcmc, burnin, df, verbose, myframe));
    return __result;
END_RCPP
}
// mlChib_cpp
double mlChib_cpp(Function lfun, arma::vec hatPar, arma::mat V, arma::mat mcmc, double df, int verbose, Environment myframe);
RcppExport SEXP iLaplace_mlChib_cpp(SEXP lfunSEXP, SEXP hatParSEXP, SEXP VSEXP, SEXP mcmcSEXP, SEXP dfSEXP, SEXP verboseSEXP, SEXP myframeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Function >::type lfun(lfunSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatPar(hatParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mcmc(mcmcSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Environment >::type myframe(myframeSEXP);
    __result = Rcpp::wrap(mlChib_cpp(lfun, hatPar, V, mcmc, df, verbose, myframe));
    return __result;
END_RCPP
}
// mlIS_cpp
double mlIS_cpp(Function lfun, arma::vec hatPar, arma::mat V, int nsim, int p, double df, int verbose, Environment myframe);
RcppExport SEXP iLaplace_mlIS_cpp(SEXP lfunSEXP, SEXP hatParSEXP, SEXP VSEXP, SEXP nsimSEXP, SEXP pSEXP, SEXP dfSEXP, SEXP verboseSEXP, SEXP myframeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Function >::type lfun(lfunSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hatPar(hatParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Environment >::type myframe(myframeSEXP);
    __result = Rcpp::wrap(mlIS_cpp(lfun, hatPar, V, nsim, p, df, verbose, myframe));
    return __result;
END_RCPP
}
// nlpost_bod2
double nlpost_bod2(arma::vec beta, double lsig, arma::vec y, arma::vec x, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_nlpost_bod2(SEXP betaSEXP, SEXP lsigSEXP, SEXP ySEXP, SEXP xSEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(nlpost_bod2(beta, lsig, y, x, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// nlpostT_bod2
double nlpostT_bod2(arma::vec beta, double lsig, double lnu, arma::vec y, arma::vec x, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_nlpostT_bod2(SEXP betaSEXP, SEXP lsigSEXP, SEXP lnuSEXP, SEXP ySEXP, SEXP xSEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< double >::type lnu(lnuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(nlpostT_bod2(beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// grad_bod2
arma::vec grad_bod2(arma::vec beta, double lsig, arma::vec y, arma::vec x, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_grad_bod2(SEXP betaSEXP, SEXP lsigSEXP, SEXP ySEXP, SEXP xSEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(grad_bod2(beta, lsig, y, x, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// hess_bod2
arma::mat hess_bod2(arma::vec beta, double lsig, arma::vec y, arma::vec x, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_hess_bod2(SEXP betaSEXP, SEXP lsigSEXP, SEXP ySEXP, SEXP xSEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(hess_bod2(beta, lsig, y, x, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// gradT_bod2
arma::vec gradT_bod2(arma::vec beta, double lsig, double lnu, arma::vec y, arma::vec x, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_gradT_bod2(SEXP betaSEXP, SEXP lsigSEXP, SEXP lnuSEXP, SEXP ySEXP, SEXP xSEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< double >::type lnu(lnuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(gradT_bod2(beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// hessT_bod2
arma::mat hessT_bod2(arma::vec beta, double lsig, double lnu, arma::vec y, arma::vec x, int n, arma::vec muBeta, arma::mat SigBeta, double sigScale);
RcppExport SEXP iLaplace_hessT_bod2(SEXP betaSEXP, SEXP lsigSEXP, SEXP lnuSEXP, SEXP ySEXP, SEXP xSEXP, SEXP nSEXP, SEXP muBetaSEXP, SEXP SigBetaSEXP, SEXP sigScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lsig(lsigSEXP);
    Rcpp::traits::input_parameter< double >::type lnu(lnuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muBeta(muBetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SigBeta(SigBetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigScale(sigScaleSEXP);
    __result = Rcpp::wrap(hessT_bod2(beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale));
    return __result;
END_RCPP
}
// seqMat
List seqMat(arma::vec par, arma::vec se, int lengthOut, int q, double delta);
RcppExport SEXP iLaplace_seqMat(SEXP parSEXP, SEXP seSEXP, SEXP lengthOutSEXP, SEXP qSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type se(seSEXP);
    Rcpp::traits::input_parameter< int >::type lengthOut(lengthOutSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    __result = Rcpp::wrap(seqMat(par, se, lengthOut, q, delta));
    return __result;
END_RCPP
}
// ldetHessBlocks
arma::vec ldetHessBlocks(arma::mat hessMat, int dimMat);
RcppExport SEXP iLaplace_ldetHessBlocks(SEXP hessMatSEXP, SEXP dimMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type hessMat(hessMatSEXP);
    Rcpp::traits::input_parameter< int >::type dimMat(dimMatSEXP);
    __result = Rcpp::wrap(ldetHessBlocks(hessMat, dimMat));
    return __result;
END_RCPP
}
// SEv
arma::vec SEv(arma::mat hessMat, int dimMat);
RcppExport SEXP iLaplace_SEv(SEXP hessMatSEXP, SEXP dimMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type hessMat(hessMatSEXP);
    Rcpp::traits::input_parameter< int >::type dimMat(dimMatSEXP);
    __result = Rcpp::wrap(SEv(hessMat, dimMat));
    return __result;
END_RCPP
}
// nlpost_rebin
double nlpost_rebin(arma::vec u, arma::vec theta, arma::vec y, int n);
RcppExport SEXP iLaplace_nlpost_rebin(SEXP uSEXP, SEXP thetaSEXP, SEXP ySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(nlpost_rebin(u, theta, y, n));
    return __result;
END_RCPP
}
// grU_rebin
arma::vec grU_rebin(arma::vec u, arma::vec theta, arma::vec y, int n);
RcppExport SEXP iLaplace_grU_rebin(SEXP uSEXP, SEXP thetaSEXP, SEXP ySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(grU_rebin(u, theta, y, n));
    return __result;
END_RCPP
}
// hessU_rebin
arma::mat hessU_rebin(arma::vec u, arma::vec theta, arma::vec y, int n);
RcppExport SEXP iLaplace_hessU_rebin(SEXP uSEXP, SEXP thetaSEXP, SEXP ySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(hessU_rebin(u, theta, y, n));
    return __result;
END_RCPP
}
// nlogH
double nlogH(arma::vec randPar, arma::vec fixPar, List data);
RcppExport SEXP iLaplace_nlogH(SEXP randParSEXP, SEXP fixParSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type randPar(randParSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fixPar(fixParSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    __result = Rcpp::wrap(nlogH(randPar, fixPar, data));
    return __result;
END_RCPP
}
// grad_nlogH
arma::vec grad_nlogH(arma::vec randPar, arma::vec fixPar, List data);
RcppExport SEXP iLaplace_grad_nlogH(SEXP randParSEXP, SEXP fixParSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type randPar(randParSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fixPar(fixParSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    __result = Rcpp::wrap(grad_nlogH(randPar, fixPar, data));
    return __result;
END_RCPP
}
// hess_nlogH
arma::mat hess_nlogH(arma::vec randPar, arma::vec fixPar, List data);
RcppExport SEXP iLaplace_hess_nlogH(SEXP randParSEXP, SEXP fixParSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type randPar(randParSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fixPar(fixParSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    __result = Rcpp::wrap(hess_nlogH(randPar, fixPar, data));
    return __result;
END_RCPP
}
// ldetHess_nlogH
NumericVector ldetHess_nlogH(arma::vec randPar, arma::vec fixPar, List data);
RcppExport SEXP iLaplace_ldetHess_nlogH(SEXP randParSEXP, SEXP fixParSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type randPar(randParSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fixPar(fixParSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    __result = Rcpp::wrap(ldetHess_nlogH(randPar, fixPar, data));
    return __result;
END_RCPP
}
