#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec SEv(arma::mat hessMat, int dimMat) {
  arma::vec out(dimMat, arma::fill::zeros);
  out(dimMat-1) = sqrt(1.0/hessMat(dimMat-1,dimMat-1));
  //  arma::mat invMat(dimMat, dimMat);
  for(int i=0; i<(dimMat-1); i++){
    arma::mat invMat(dimMat-i, dimMat-i);
    invMat = inv_sympd(hessMat(arma::span(i,dimMat-1),arma::span(i,dimMat-1)));
    out(i) = sqrt(invMat(0,0));
  }
  return out;
}

// [[Rcpp::export]]
arma::vec seq_C(double a, double cent, double b, int lengthOut) {
  arma::vec out(lengthOut);
  out.fill(a);
  out[lengthOut-1] = b;
  double width = (b-a)/(lengthOut-1);
  for(int i=1; i<(lengthOut -1); i++){
    out(i) = out(i-1) + width;
  }
  out(lengthOut/2) = cent;
  return out;
}

// [[Rcpp::export]]
List seqMat(arma::vec par, arma::vec se, int lengthOut, int q, double delta){
  arma::vec av = par - delta*se;
  arma::vec bv = par + delta*se;
  arma::mat out(lengthOut+1, q, arma::fill::zeros);
  for(int i=0; i<q; i++){
    out(arma::span(0,lengthOut-1),i) = seq_C(av(i), par(i), bv(i), lengthOut);
    out(lengthOut,i) = i+1;
  }

  return List::create(Named("parVal") = out,
                      Named("lo") = av,
                      Named("up") = bv);
}

// [[Rcpp::export]]
arma::vec ldetHessBlocks(arma::mat hessMat, int dimMat){
  arma::vec out(dimMat-1, arma::fill::zeros);
  double ldet = 0.0, sign = 1.0;
  for(int i=0; i<(dimMat-1); i++){
    arma::log_det(ldet, sign, hessMat(arma::span(i,dimMat-1),
                                      arma::span(i,dimMat-1)));
    out(i) = ldet;
  }
  return out;
}
