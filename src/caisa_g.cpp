#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// FUNCTION DECLARATIONS
// ---------------------
arma::mat outerAddition_g     (const arma::vec& a, const arma::vec& b);
void updatebetaj_g            (const arma::vec& xj, double wj,
                               double& betaj, arma::vec& r,
                               arma::mat& phi, arma::vec& pi, arma::vec& piold,
                               double sigma2, const arma::vec& sa2,
                               const arma::vec& s2inv,
                               double& a1, double& a2,
                               int j, int p, bool flag,
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
                       double stepsize, bool updatesigma, bool mode,
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
  arma::mat S2inv        = 1 / outerAddition_g(1/sa2, w);
  S2inv.row(0).fill(epstol);
  
  // ---------------------------------------------------------------------
  // INITIALIZE
  // ---------------------------------------------------------------------
  int iter               = 0;
  int j;
  arma::vec varobj(maxiter);
  double a1;
  double a2;
  bool flag              = true;
  arma::vec piold;
  
  // ---------------------------------------------------------------------
  // START LOOP : CYCLE THROUGH COORDINATE ASCENT UPDATES
  // ---------------------------------------------------------------------
  for (iter = 0; iter < maxiter; iter++) {
    
    // reset parameters
    a1                   = 0;
    a2                   = 0;
    piold                = pi;
    if (!flag) {
      pi.fill(0);
    }
    
    // ---------------------------------------------------------------------
    // RUN COORDINATE ASCENT UPDATES : INDEX 1 - INDEX P
    // ---------------------------------------------------------------------
    for (j = 0; j < p; j++){
      
      updatebetaj_g(X.col(j), w(j), beta(j), r, phi, pi, piold, sigma2, sa2,
		                S2inv.col(j), a1, a2, j, p, flag, stepsize, epstol);
      
    }
    
    // ---------------------------------------------------------------------
    // MAKE SURE PI IS POSITIVE
    // ---------------------------------------------------------------------
    if (flag) {
      pi.elem(find(pi < 0)).fill(0);
      pi                  = pi / sum(pi);
    }
    
    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE
    // ---------------------------------------------------------------------
    varobj(iter)          = dot(r,r)/2 - dot(square(beta), w)/2 + a1/2;
    
    if (updatesigma)
      sigma2              = 2 * varobj(iter) / n;
    
    varobj(iter)          = varobj(iter) / sigma2 + log(2*PI*sigma2)/2 * n -
                            dot(pi, log(pi + epstol)) * p + a2;
    
    for (j = 1; j < K; j++){
      varobj(iter)       +=  pi(j) * log(sa2(j)) * p / 2;
    }
    
    // ---------------------------------------------------------------------
    // CHECK CONVERGENCE
    // ---------------------------------------------------------------------
    if (iter >= miniter - 1) {
      if (max(abs(pi - piold)) < convtol * K) {
        if (mode) {
          iter++;
          break;
        } else {
          flag            = false;
          mode            = true;
        }
      }
      
      if (iter > 0) {
        if (varobj(iter) > varobj(iter - 1)){
          if (mode) {
            iter++;
            break;
          } else {
            flag           = false;
            mode           = true;
          }
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

void updatebetaj_g     (const arma::vec& xj, double wj,
                        double& betaj, arma::vec& r,
                        arma::mat& phi, arma::vec& pi, arma::vec& piold,
                        double sigma2, const arma::vec& sa2,
                        const arma::vec& s2inv,
                        double& a1, double& a2,
                        int j, int p, bool flag,
                        double stepsize, double epstol) {
  
  // update pi if flag
  if (flag) {
    pi                 += -trans(phi.row(j)) / p;
  }
  
  // calculate b
  double bjwj           = dot(r, xj) + betaj * wj;
  
  // update r first step
  r                    += xj * betaj;
  
  // calculate muj
  arma::vec muj         = bjwj * s2inv;
  arma::vec phij;
  
  if (flag) {
    // calculate Lj
    arma::vec Lj        = -log(1 + sa2 * wj)/2 + muj * (bjwj / 2 / sigma2);
    Lj                  = exp(Lj - max(Lj));
    
    // calculate phij
    phij                = (pi + epstol) % Lj;
    phij               %= arma::square(.5 + stepsize * (Lj / p + dot(Lj, pi)));
    phij                = phij / sum(phij);
    phi.row(j)          = trans(phij);
    
  } else {
    // calculate phij
    phij                = log(piold + epstol) - log(1 + sa2 * wj)/2 + muj * (bjwj / 2 / sigma2);
    phij                = exp(phij - max(phij));
    phij                = phij / sum(phij);
  }
  
  // update pi
  pi                   += phij / p;
  
  // update betaj
  betaj                 = dot(phij, muj);
  
  // update r second step
  r                    += -xj * betaj;
  
  // precalculate for M-step
  a1                   += bjwj * betaj;
  a2                   += dot(phij, log(phij + epstol));
  phij(0)               = 0;
  a2                   += -dot(phij, log(s2inv)) / 2;
  
  return;
}

arma::mat outerAddition_g (const arma::vec& a, const arma::vec& b) {
  arma::mat A(a.n_elem, b.n_elem);
  A.fill(0);
  A.each_row() += b.t();
  A.each_col() += a;
  return A;
}
