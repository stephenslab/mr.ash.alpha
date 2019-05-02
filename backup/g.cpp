#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// FUNCTION DECLARATIONS
// ---------------------
arma::mat outerAddition     (const arma::vec& a, const arma::vec& b);
void updatebetaj            (const arma::vec& xj, double wj,
                             double& betaj, arma::vec& r,
                             arma::mat& phi, arma::vec& pi,
                             double sigma2, const arma::vec& sa2,
                             const arma::vec& s2inv,
                             arma::vec& a1, arma::vec& a2,
                             int j, int p,
                             double stepsize, double epstol);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
List caisa_g          (const arma::mat& X,
                       const arma::vec& w, const arma::vec& sa2,
                       arma::mat& phi,
                       arma::vec& pi, arma::vec& beta,
                       arma::vec& r, double sigma2,
                       int maxiter, int miniter, double convtol, double epstol,
                       double stepsize, bool updatesigma,
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
  arma::mat S2inv        = 1 / outerAddition(1/sa2, w);
  S2inv.row(0).fill(epstol);
  
  
  // ---------------------------------------------------------------------
  // INITIALIZE
  // ---------------------------------------------------------------------
  int iter               = 0;
  int j;
  arma::vec varobj(maxiter);
  arma::vec a1(p);
  arma::vec a2(p);
  arma::vec piold;
  
  // ---------------------------------------------------------------------
  // START LOOP : CYCLE THROUGH COORDINATE ASCENT UPDATES
  // ---------------------------------------------------------------------
  for (iter = 0; iter < maxiter; iter++) {
    
    piold                = pi;
    
    // ---------------------------------------------------------------------
    // RUN COORDINATE ASCENT UPDATES : INDEX 1 - INDEX P
    // ---------------------------------------------------------------------
    for (j = 0; j < p; j++){
      
      updatebetaj(X.col(j), w(j), beta(j), r, phi, pi, sigma2, sa2,
                  S2inv.col(j), a1, a2, j, p, stepsize, epstol);
      
    }
    
    // ---------------------------------------------------------------------
    // MAKE SURE PI IS POSITIVE
    // ---------------------------------------------------------------------
    pi.elem(find(pi < 0)).fill(0);
    pi                    = pi / sum(pi);
    
    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE
    // ---------------------------------------------------------------------
    varobj(iter)          = dot(r,r)/2 - dot(square(beta), w)/2 + sum(a1)/2;
    
    if (updatesigma)
      sigma2              = 2 * varobj(iter) / n;
    
    varobj(iter)          = varobj(iter) / sigma2 + log(2*PI*sigma2)/2 * n -
                            dot(pi, log(pi + epstol)) * p + sum(a2);
    
    for (j = 1; j < K; j++){
      varobj(iter)       +=  pi(j) * log(sa2(j)) * p / 2;
    }
    
    // ---------------------------------------------------------------------
    // CHECK CONVERGENCE
    // ---------------------------------------------------------------------
    if (iter >= miniter - 1) {
      if (max(abs(pi - piold)) < convtol * K) {
        iter++;
        break;
      }
      
      if (iter > 0) {
        if (varobj(iter) > varobj(iter - 1)){
          break;
        }
      }
    }
  }
  
  // ---------------------------------------------------------------------
  // RETURN VALUES
  // ---------------------------------------------------------------------
  
  return List::create(Named("beta")    = beta,
                      Named("sigma2")  = sigma2,
                      Named("pi")      = pi,
                      Named("iter")    = iter,
                      Named("varobj")  = varobj.subvec(0,iter-1));
}

void updatebetaj       (const arma::vec& xj, double wj,
                        double& betaj, arma::vec& r,
                        arma::mat& phi, arma::vec& pi,
                        double sigma2, const arma::vec& sa2,
                        const arma::vec& s2inv,
                        arma::vec& a1, arma::vec& a2,
                        int j, int p,
                        double stepsize, double epstol) {
  
  // update pi
  pi                   += -trans(phi.row(j)) / p;
  
  // calculate b
  double bjwj           = dot(r, xj) + betaj * wj;
  
  // update r first step
  r                    += xj * betaj;
  
  // calculate muj
  arma::vec muj         = bjwj * s2inv;
  
  // calculate Lj
  arma::vec Lj          = -log(1 + sa2 * wj)/2 + muj * (bjwj / 2 / sigma2);
  Lj                    = exp(Lj - max(Lj));
  
  // calculate phij
  arma::vec phij        = pi % Lj;
  phij                 %= pow(.5 + stepsize * (Lj / p + dot(Lj, pi)), 2);
  phij                  = phij / sum(phij);
  phi.row(j)            = trans(phij);
  
  // update pi second step
  pi                   += phij / p;
  
  // update betaj
  betaj                 = dot(phij, muj);
  
  // update r second step
  r                    += -xj * betaj;
  
  // precalculate for M-step
  a1(j)                 = dot(phij % square(muj), 1 / s2inv);
  a2(j)                 = dot(phij, log(phij + epstol));
  phij(0)               = 0;
  a2(j)                += -dot(phij, log(s2inv)) / 2;
  
  return;
}

arma::mat outerAddition    (const arma::vec& a, const arma::vec& b) {
  arma::mat A(a.n_elem, b.n_elem);
  A.fill(0);
  A.each_row() += b.t();
  A.each_col() += a;
  return A;
}
