// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// FUNCTION DECLARATIONS
// ---------------------
arma::mat outerAddition_gs  (const arma::vec& a, const arma::vec& b);
void samplebetaj            (const arma::vec& xj, double wj,
                             double& betaj, arma::vec& r,
                             arma::vec& pi,
                             double sigma2, const arma::vec& sa2,
                             const arma::vec& s2inv, int K);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
List gibbs_sampling   (const arma::mat& X,
                       const arma::vec& w, const arma::vec& sa2,
                       arma::vec& pi, arma::vec& beta,
                       arma::vec& r, double sigma2,
                       int maxiter, int burnin,
                       bool verbose) {
  
  // ---------------------------------------------------------------------
  // DEFINE SIZES
  // ---------------------------------------------------------------------
  int n                  = X.n_rows;
  int p                  = X.n_cols;
  int K                  = sa2.n_elem;
  
  // ---------------------------------------------------------------------
  // PRECALCULATE
  // ---------------------------------------------------------------------
  arma::mat S2inv        = 1 / outerAddition_gs(1/sa2, w);
  arma::vec beta_save(p, arma::fill::zeros);
  double sigma2_save = 0;
  
  // ---------------------------------------------------------------------
  // INITIALIZE
  // ---------------------------------------------------------------------
  int iter               = 0;
  int j;
  
  // ---------------------------------------------------------------------
  // START LOOP : CYCLE THROUGH GIBBS SAMPLING
  // ---------------------------------------------------------------------
  for (iter = 0; iter < maxiter; iter++) {
    
    // ---------------------------------------------------------------------
    // SAMPLE BETA_j
    // ---------------------------------------------------------------------
    for (j = 0; j < p; j++){
      
      samplebetaj(X.col(j), w(j), beta(j), r, pi, sigma2, sa2,
                  S2inv.col(j), K);
      
    }
    
    // ---------------------------------------------------------------------
    // SAMPLE SIGMA2
    // ---------------------------------------------------------------------
    sigma2              = 1 / R::rgamma( n / 2.0, 2.0 / arma::dot(r,r) );
    
    // ---------------------------------------------------------------------
    // SAVE
    // ---------------------------------------------------------------------
    if (iter >= burnin) {
      beta_save        += beta;
      sigma2_save      += sigma2;
    }
  }
  
  // ---------------------------------------------------------------------
  // RETURN VALUES
  // ---------------------------------------------------------------------
  return List::create(Named("beta")    = beta_save / (maxiter - burnin),
                      Named("sigma2")  = sigma2_save / (maxiter - burnin));
}

void samplebetaj       (const arma::vec& xj, double wj,
                        double& betaj, arma::vec& r,
                        arma::vec& pi,
                        double sigma2, const arma::vec& sa2,
                        const arma::vec& s2inv, int K) {
  
  // calculate b
  double bjwj           = dot(r, xj) + betaj * wj;
  
  // update r first step
  r                    += xj * betaj; 
  
  // calculate muj
  arma::vec muj         = bjwj * s2inv;
  
  // calculate phij
  arma::vec phij        = -log(1 + sa2 * wj)/2.0 + muj * (bjwj / 2.0 / sigma2);
  phij                  = pi % exp(phij - max(phij));
  phij                  = phij / sum(phij);
  
  // sample muj
  muj                  += sqrt(sigma2) * arma::sqrt(s2inv) % arma::randn(K);
  
  // update betaj
  betaj                 = dot(phij, muj);
  
  // update r second step
  r                    += -xj * betaj;
  
  return;
}

arma::mat outerAddition_gs (const arma::vec& a, const arma::vec& b) {
  arma::mat A(a.n_elem, b.n_elem);
  A.fill(0);
  A.each_row() += b.t();
  A.each_col() += a;
  return A;
}
