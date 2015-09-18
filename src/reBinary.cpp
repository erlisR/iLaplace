#include <RcppArmadillo.h>
using namespace Rcpp;

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

// [[Rcpp::export]]
arma::vec SEv(arma::mat hessMat, int dimMat) {
  arma::vec out(dimMat, arma::fill::zeros);
  out(dimMat-1) = sqrt(1.0/hessMat(dimMat-1,dimMat-1));
  //  arma::mat invMat(dimMat, dimMat);
  for(int i=0; i<(dimMat-1); i++){
    arma::mat invMat(dimMat-i, dimMat-i);
    invMat = inv(hessMat(arma::span(i,dimMat-1),arma::span(i,dimMat-1)));
    out(i) =   sqrt(invMat(0,0));
  }
  return out;
}

// [[Rcpp::export]]
double nlpost_rebin (arma::vec u, arma::vec theta, arma::vec y, int n){
  double eta, pbin, ll, pr;
  ll = 0.0;
  for(int i=0; i<n; i++){
    eta = theta(0) + u(i);
    pbin = 1.0/(1.0 + exp(-eta));
    // lpbin = -log(exp(0.0) + exp(-pbin));
//     
//     ll += y(i)*lpbin + (1.0-y(i))*log(1.0-pbin) +
//         Rf_dnorm4(u(i), 0.0, 1.0/sqrt(exp(theta(1))),true);
    ll += Rf_dbinom(y(i),1, pbin, true) + 
          Rf_dnorm4(u(i), 0.0, 1.0/sqrt(exp(theta(1))),true);
  }
  // Rf_dgamma(exp(theta(1)), 1.0, 1.0, true) + theta(1)
  pr = Rf_dgamma(exp(theta(1)), 1.0, 1.0, true) +
       theta(1) +
       Rf_dnorm4(theta(0),0.0, 1.0, true);
  return -(ll + pr);
}

// [[Rcpp::export]]
arma::vec grU_rebin (arma::vec u, arma::vec theta, arma::vec y, int n){
  arma::vec ans(n);
  double eta=0.0, etaPrime=0.0;
  
  for(int i=0; i<n; i++){
    eta = 1.0/(1.0+exp(-1.0*(theta(0)+u(i))));
    etaPrime = eta*eta*exp(-1.0*(theta(0)+u(i)));
    
    ans(i) = y(i)*etaPrime/eta - 
             (1-y(i))*etaPrime/(1-eta) -
             u(i)*exp(theta(1));
  }
  return -ans;
}


arma::vec grUB_rebin(arma::vec u, arma::vec theta, arma::vec y, int n){
  arma::vec ans(n+1);
  double eta  = 0.0, etaPrimeu = 0.0,
    grbeta = 0.0;
  for(int i=0; i<n; i++){
    eta = 1.0/(1.0+exp(-1.0*(theta(0) + u(i))));
    etaPrimeu = eta*eta*exp(-1.0*(theta(0)+u(i)));
    ans(i) = y(i)*etaPrimeu/eta -
      (1-y(i))*etaPrimeu/(1-eta) -
      u(i)*exp(theta(1));
    grbeta += y(i)*etaPrimeu/eta - (1-y(i))*etaPrimeu/(1- eta);
  }
  ans(n) = grbeta -theta(0);
  // out(0) = ans(n);
  // out(arma::span(1,n)) = ans(arma::span(0,(n-1)));
  return -ans;
}


arma::vec grAll_rebin (arma::vec u, arma::vec theta, arma::vec y, int n){
  arma::vec ans(n), out(n+2);
  double eta  = 0.0, etaPrimeu = 0.0,
    grbeta = 0.0, grlsig2 = 0.0;
  for(int i=0; i<n; i++){
    eta = 1.0/(1.0+exp(-1.0*(theta(0) + u(i))));
    etaPrimeu = eta*eta*exp(-1.0*(theta(0)+u(i)));
    ans(i) = y(i)*etaPrimeu/eta -
              (1-y(i))*etaPrimeu/(1-eta) -
              u(i)*exp(theta(1));
    grbeta += y(i)*etaPrimeu/eta - (1-y(i))*etaPrimeu/(1- eta);
    //ricorda di aggiungere alla fine la derivata della prior di beta!!
    grlsig2 += 0.5*(1.0 - u(i)*u(i)*exp(theta(1)));
      //ricorda di aggiungere alla fine la derivata della prior di beta!!
  }
  // (u, b, lprec)
//   out(n+1) = grlsig2 - exp(theta(1)) + 1.0;
//   out(n) = grbeta -theta(0);
//   out(arma::span(0,n-1))=ans;

  out(1) = grlsig2 - exp(theta(1)) + 1.0;
  out(0) = grbeta -theta(0);
  out(arma::span(2,n+1))=ans;
  
  return -out;
}

// [[Rcpp::export]]
arma::mat hessU_rebin (arma::vec u, arma::vec theta, arma::vec y, int n){
  arma::mat ans(n, n, arma::fill::zeros);
  double eta=0.0, etaPrime=0.0, etaSecnd=0.0;
  for(int i=0; i<n; i++){
    eta = 1.0/(1.0+exp(-1.0*(theta(0) + u(i))));
    etaPrime = eta*eta*exp(-1.0*(theta(0)+u(i)));
    etaSecnd = 2*exp(-2*(theta(0)+u(i)))*pow(eta,3.0)-etaPrime;
    
    ans(i,i) = -(1-y(i))*pow(etaPrime/(1-eta),2.0) -
                y(i)*pow(etaPrime/eta,2.0) -
                (1-y(i))*etaSecnd/(1-eta) +
                y(i)*etaSecnd/eta - 
                exp(theta(1));
  }
  return -ans;
}

arma::mat hessUB_rebin(arma::vec u, arma::vec theta, arma::vec y, int n){
  arma::mat uu(n, n, arma::fill::zeros), ans(n+1,n+1,arma::fill::zeros);
  arma::vec bu(n);
  arma::mat out(n+1,n+1);
  double bb = 0.0, eta  = 0.0, etaPrimeu = 0.0, etaSecndu=0.0;
  for(int i=0; i<n; i++){
    eta = 1.0/(1.0+exp(-1.0*(theta(0) + u(i))));
    etaPrimeu = eta*eta*exp(-1.0*(theta(0)+u(i)));
    etaSecndu = 2*exp(-2*(theta(0)+u(i)))*pow(eta,3.0)-etaPrimeu;
    uu(i,i) = -(1.0 + exp(2.0*(theta(0) + u(i))) + 
      exp(theta(0) + u(i))*(2.0 + 1.0/exp(theta(1))))*exp(theta(1))/
    pow(1.0 + exp(theta(0) + u(i)), 2.0);
    bb += -(1-y(i))*pow(etaPrimeu/(1-eta),2.0) -
      y(i)*pow(etaPrimeu/eta,2.0) -
      (1-y(i))*etaSecndu/(1-eta) +
      y(i)*etaSecndu/eta;
    //ricorda di aggiungere a bb la derivata seconda della logprior
    bu(i) = y(i)*(etaSecndu*eta-etaPrimeu*etaPrimeu)/pow(eta, 2.0) -
      (1-y(i))*(etaSecndu*(1-eta) + etaPrimeu*etaPrimeu)/pow(1-eta,2.0);
  }
  ans(arma::span(0,n-1),arma::span(0,n-1)) = uu;
  ans(n,n) = bb - 1;
  ans(n, arma::span(0,n-1)) = bu.t();
  ans(arma::span(0,n-1), n) = bu;
  //   out(0,0) = bb - 1;
  //   out(0,arma::span(1,n)) = bu.t();
  //   out(arma::span(1,n),arma::span(1,n)) = uu;
  //   out(arma::span(1,n), 0) = bu;
  return -ans;
}


arma::mat hessAll_rebin (arma::vec u, arma::vec theta, arma::vec y, int n){
  arma::mat uu(n, n, arma::fill::zeros), ans(n+2,n+2,arma::fill::zeros);
  arma::vec bu(n), su(n);
  double bb = 0.0, ss = 0.0, eta  = 0.0, etaPrimeu = 0.0, etaSecndu=0.0;
  for(int i=0; i<n; i++){
    eta = 1.0/(1.0+exp(-1.0*(theta(0) + u(i))));
    etaPrimeu = eta*eta*exp(-1.0*(theta(0)+u(i)));
    etaSecndu = 2*exp(-2*(theta(0)+u(i)))*pow(eta,3.0)-etaPrimeu;
    
    uu(i,i) = -(1-y(i))*pow(etaPrimeu/(1-eta),2.0) -
                y(i)*pow(etaPrimeu/eta,2.0) -
                (1-y(i))*etaSecndu/(1-eta) +
                y(i)*etaSecndu/eta - 
                exp(theta(1));
    
    bb += -(1-y(i))*pow(etaPrimeu/(1-eta),2.0) -
           y(i)*pow(etaPrimeu/eta,2.0) -
           (1-y(i))*etaSecndu/(1-eta) +
            y(i)*etaSecndu/eta;
    
    //ricorda di aggiungere a bb la derivata seconda della logprior
    ss -= 0.5*u(i)*u(i)*exp(theta(1));
    bu(i) = y(i)*(etaSecndu*eta-etaPrimeu*etaPrimeu)/pow(eta, 2.0) -
            (1-y(i))*(etaSecndu*(1-eta) + etaPrimeu*etaPrimeu)/pow(1-eta,2.0);
    //     bu(i) = -(1-y(i))*etaPrimeu*etaPrimeu/pow(1-eta,2.0)-
//             y(i)*etaPrimeu*etaPrimeu/(eta*eta)-
//             (1-y(i))*etaSecndu/(1-eta)+
//             y(i)*etaSecndu/eta;
    su(i) = -u(i)*exp(theta(1));
  }
// ordering (b, tau, u)
  ans(arma::span(2,n+1),arma::span(2,n+1)) = uu;
  ans(0,0) = bb - 1;
  ans(1,1) = ss - exp(theta(1));
  ans(0,1) = ans(1,0) = 0.0;
  ans(0, arma::span(2,n+1)) = bu.t();
  ans(1, arma::span(2,n+1)) = su.t();
  ans(arma::span(2,n+1), 0) = bu;
  ans(arma::span(2,n+1), 1) = su;

//   ordering (u, b, lprec)
//   ans(arma::span(0,n-1),arma::span(0,n-1)) = uu;
//   ans(n,n) = bb - 1;
//   ans(n+1,n+1) = ss - exp(theta(1));
//   ans(n,n+1) = ans(n+1,n) = 0.0;
//   ans(n, arma::span(0,n-1)) = bu.t();
//   ans(n+1, arma::span(0,n-1)) = su.t();
//   ans(arma::span(0,n-1), n) = bu;
//   ans(arma::span(0,n-1), n+1) = su;
  return -ans;
}
