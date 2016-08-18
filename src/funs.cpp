#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List aux_quant(arma::mat V, int nc) {

  // out is a list with nc+1 elements, where each
  // element is a vector of various dimensions. In
  // particular:
  // - out(0) contains the vector of std.dev of of the various densities.
  // - out(1) to out(nc-1) contain the vectors of the regression
  // coefs. useful for the nested optimisations.
  // - out(nc) contains the vector of log-determinants of  blocks of V

  Rcpp::List out(nc+1); // the elements of the list are left unspecified
  arma::mat S = arma::inv_sympd(V);
  arma::vec se(nc, arma::fill::zeros), ldetb(nc, arma::fill::zeros);

  // int j = 0;
  double sign = 1.0;
  se(0) = S(0,0);
  arma::log_det(ldetb(0), sign, V);
  ldetb(nc-1) = log(V(nc-1,nc-1));

  for(int i=1; i<(nc-1); i++){
    se(i) = arma::as_scalar(S(i,i) -
      S(i,arma::span(0,i-1))*
      inv_sympd(S(arma::span(0,i-1),arma::span(0,i-1)))*
      trans(S(i,arma::span(0,i-1))));
    arma::log_det(ldetb(i), sign, V(arma::span(i,nc-1),arma::span(i,nc-1)));

    // out(i) = arma::as_scalar(V(arma::span(i)))
    // cnst(j*,j*);
    out(i) = inv(V(arma::span(i,nc-1),arma::span(i,nc-1)))*V(arma::span(i,nc-1),i-1);
    // std::cout << cnst << std::endl;

  }
  se(nc-1) = 1.0/V(nc-1,nc-1);
  out(0) = sqrt(se);
  out(nc-1) = V(nc-1,nc-2)/V(nc-1,nc-1);
  out(nc) = ldetb;
  // out.names() = "se";
  // out.names() = "ldetb";
  return out;
}

// // [[Rcpp::export]]
// arma::vec SEv(arma::mat hessMat, int dimMat) {
//   arma::vec out(dimMat, arma::fill::zeros);
//   out(dimMat-1) = sqrt(1.0/hessMat(dimMat-1,dimMat-1));
//   //  arma::mat invMat(dimMat, dimMat);
//   for(int i=0; i<(dimMat-1); i++){
//     arma::mat invMat(dimMat-i, dimMat-i);
//     invMat = inv_sympd(hessMat(arma::span(i,dimMat-1),arma::span(i,dimMat-1)));
//     out(i) = sqrt(invMat(0,0));
//   }
//   return out;
// }
//
// // [[Rcpp::export]]
// Rcpp::List fnblocks(arma::mat hessMat, int dimMat) {
//   arma::vec sev(dimMat, arma::fill::zeros);
//   arma::vec ldetv(dimMat, arma::fill::zeros);
//   Rcpp::List out(dimMat+2);
//
//   double sign = 1.0;
//
//   //  arma::mat invMat(dimMat, dimMat);
//   out(0) = 0.0;
//   for(int i=0; i<(dimMat-1); i++){
//     arma::mat invMat(dimMat-i, dimMat-i);
//     invMat = inv_sympd(hessMat(arma::span(i,dimMat-1),arma::span(i,dimMat-1)));
//     sev(i) = sqrt(invMat(0,0));
//     if(i>0) {
//       out(i) = invMat*hessMat(arma::span(i,dimMat-1),i-1);
//     }
//     arma::log_det(ldetv(i), sign,
//                   hessMat(arma::span(i,dimMat-1),arma::span(i,dimMat-1)));
//   }
//   ldetv(dimMat-1) = hessMat(dimMat-1,dimMat-1);
//   sev(dimMat-1) = sqrt(1.0/hessMat(dimMat-1,dimMat-1));
//   out(dimMat-1) = hessMat(dimMat-1,dimMat-2)/hessMat(dimMat-1,dimMat-1);
//   out(dimMat) = sev;
//   out(dimMat+1) = ldetv;
//   return out;
// }
//
// // [[Rcpp::export]]
// arma::vec seq_C(double a, double cent, double b, int lengthOut) {
//   arma::vec out(lengthOut);
//   out.fill(a);
//   out[lengthOut-1] = b;
//   double width = (b-a)/(lengthOut-1);
//   for(int i=1; i<(lengthOut -1); i++){
//     out(i) = out(i-1) + width;
//   }
//   out(lengthOut/2) = cent;
//   return out;
// }
//
// // [[Rcpp::export]]
// Rcpp::List seqMat(arma::vec par, arma::vec se, int lengthOut, int q, double delta){
//   arma::vec av = par - delta*se;
//   arma::vec bv = par + delta*se;
//   arma::mat out(lengthOut+1, q, arma::fill::zeros);
//   for(int i=0; i<q; i++){
//     out(arma::span(0,lengthOut-1),i) = seq_C(av(i), par(i), bv(i), lengthOut);
//     out(lengthOut,i) = i+1;
//   }
//
//   return Rcpp::List::create(Rcpp::Named("parVal") = out,
//                             Rcpp::Named("lo") = av,
//                             Rcpp::Named("up") = bv);
// }
//
// // [[Rcpp::export]]
// arma::vec ldetHessBlocks(arma::mat hessMat, int dimMat){
//   arma::vec out(dimMat-1, arma::fill::zeros);
//   double ldet = 0.0, sign = 1.0;
//   for(int i=0; i<(dimMat-1); i++){
//     arma::log_det(ldet, sign, hessMat(arma::span(i,dimMat-1),
//                                       arma::span(i,dimMat-1)));
//     out(i) = ldet;
//   }
//   return out;
// }
