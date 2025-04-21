#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec p_correction(const arma::colvec& pa,
                          const arma::colvec& p0,
                          const bool& log10p = true){
  
  arma::uvec s = arma::sort_index(pa);
  arma::colvec pa_sort = pa(s);
  s = arma::sort_index(s);
  arma::colvec p0_sort = arma::sort(p0);
  arma::colvec pc(pa.n_rows, arma::fill::ones);
  
  int j = 0;
  for(int i = 0; i < pc.n_rows; i++){
    
    while((j < p0.n_rows) && (p0_sort(j) < pa_sort(i))){
      
      j++;
      
    }
    pc(i) = j;
    
  }
  
  pc = pc(s);
  pc = pc / p0.n_rows;
  
  // Put into log10-scale if log10 = true:
  arma::uvec w_smaller;
  if(log10p){
    
    pc = arma::log10(pc);
    // Find where the raw p-values are smaller than 1 / (# of permuted p-values):
    w_smaller = arma::find(pa < log10(1.0 / p0.n_rows));
    
  } else {
    
    w_smaller = arma::find(pa < (1.0 / p0.n_rows));
    
  }
  
  // At those locations, we replace the corrected p-values with the raw:
  pc(w_smaller) = pa(w_smaller);
  
  return(pc);
  
}
