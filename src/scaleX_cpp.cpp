#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat scaleX_cpp(const arma::mat& X,
                     const double& tol = 1e-10){
  
  arma::mat Xc = X;
  
  // Column-center:
  Xc.each_row() -= mean(Xc, 0);
  
  // Column-scale:
  arma::rowvec s = sqrt(mean(Xc % Xc, 0));
  s.transform( [tol](double val) { return (val > tol) ? val : std::numeric_limits<double>::infinity(); } );
  
  Xc.each_row() %= (1 / s);
  
  // Return:
  return(Xc);
  
}