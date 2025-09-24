#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double sumC(const arma::mat& M,
            const arma::mat& E){
  
  int r = M.n_cols;
  int m = E.n_cols;
  
  arma::mat tE = E.t();
  arma::mat E2 = E;
  arma::mat C(m, m);
  
  double sumc = 0;
  for(int i = 0; i < r; i++){
    for(int j = 0; j <= i; j++){
      E2 = E;
      E2.each_col() %= (M.col(i) % M.col(j));
      C = tE * E2;
      if(i == j){
        sumc += arma::accu(C % C);
      } else {
        sumc += 2 * arma::accu(C % C);
      }
    }
  }
  
  return(sumc);
  
}