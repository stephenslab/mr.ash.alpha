#' @title Multiple Regression with Adaptive Shrinkage
#' 
#' @description Implements the Variational Empirical Bayes (VEB)
#' method for multiple linear regression. It maximizes the approximate
#' marginal likelihood ("evidence lower bound", or "ELBO") using a
#' coordinate ascent algorithm.
#' 
#' @details The VEB approach is based on the multiple linear
#' regression model: \deqn{y|X,\beta,\sigma^2 ~ N(X\beta, \sigma^2
#' I_n), \beta | \pi, \sigma^2 ~ \sum_{k=1}^K N(0,\sigma^2\sigma_k^2)}
#' Here \eqn{\sigma_k^2} is the k-th mixture component variance
#' \code{sa2[k]} and \eqn{K} is the number of mixture components
#' \code{length(sa2)}.  The other parameters are described in the
#' \sQuote{Arguments}.
#' 
#' The VEB approach solves the following optimization problem:
#' \deqn{F(q,g,\sigma^2) = E_q \log p(y|X,\beta,\sigma^2) -
#' \sum_{j=1}^p D_{KL}(q_j || g)} The algorithm updates the
#' variational factors \eqn{q_1,...,q_p}, \eqn{g} and \eqn{\sigma^2}
#' one at a time while fixing the others, in each outer loop
#' iteration.
#' 
#' The algorithm does not store the full variational posterior \eqn{q
#' = (q_1,...,q_p)}, but only stores the variational posterior mean
#' \code{beta} for each regression coefficients.  In order to recover
#' the full posterior, see the documnetation for
#' \code{get.full.posterior} function.
#' 
#' See \sQuote{References} for more details about the VEB approach.
#' 
#' @seealso \code{\link{get.full.posterior}}.
#'
#' @param X The input matrix, of dimension (n,p); each column is a
#' single predictor; and each row is an observation vector. Here n is
#' the number of samples and p is the number of predictors. Currently
#' sparse matrix formats are not supported.
#' 
#' @param y The response variable. Currently we only allow the linear regression case
#' which corresponds to family = "gaussian" in glmnet package.
#' Thus we treat y as a real valued quantitative response variable.
#' 
#' @param Z The covariate matrix, of dimension (n,k); k is the number of covariates.
#' The input matrix Z can be modified according to "intercept" argument.
#' If \code{Z = NULL} and \code{intercept = TRUE}, then the actual \code{Z} will be
#' the matrix having entries 1 of dimension (n,1).
#' (\code{Z = NULL} by default, and if \code{intercept = FALSE}
#' then we do not include any covariates in the model.
#' If \code{intercept = TRUE}, then we will add the vector of ones to the columns of Z.
#' That is, \code{Z <- cbind(1,Z)}.)
#' 
#' @param sa2 The vector of mixture component variances. Currently we only allow `sa2[1] = 0`.
#' for a technical reason. The default value is \code{sa2[k] = 2^(k-1) - 1}, for k = 1,...,20.
#' 
#' For Dev: (1) accelerate and block use different updates for g
#' (2) caisa, sigma, sigma_scaled, sigma_indep use different updates for sigma2,
#' based on different parametrizations
#' depending on whether sigma-dependent q, g and beta are used or not.
#' See reference for details.
#' 
#' @param method In the manuscript (preprint) listed in \sQuote{References}, only
#' \code{method = "caisa"} is used, which stands for Cooridinate Ascent Iterative Shinkage
#' Algorithm. Other method arguments will work, and produce similar outcomes unless the
#' regression setting is extreme.
#' 
#' @param max.iter The maximum number of outer loop iterations allowed.
#' 
#' @param min.iter The minimum number of inner loop iterations allowed.
#' 
#' @param beta.init The initial value for the variational posterior mean
#' of the regression coefficients.
#' 
#' @param update.pi The boolean parameter indicating whether the mixture proportion
#' \eqn{\pi} will be updated or not. In the manuscript, \code{update.pi = TRUE}.
#' 
#' @param pi The initial value for the mixture proportions \eqn{\pi_1,...,\pi_K}.
#' If \code{pi = NULL}, the default value \code{pi[k] = 1/K} for k = 1,...,K will be used.
#' 
#' @param update.sigma2 The boolean parameter indicating whether the noise variance
#' \eqn{\sigma^2} will be updated or not. In the manuscript, \code{update.sigma = TRUE}.
#' 
#' @param sigma2 The initial value for the noise variance \eqn{\sigma^2}.
#' If \code{sigma2 = NULL}, the default value \code{var(y-X\%*\%beta)} will be used.
#' 
#' @param standardize The logical flag for standardization of the columns of X variable,
#' prior to the model fitting. The coefficients are always returned on the original scale.
#' 
#' @param intercept The logical flag for including intercept (\code{intercept = TRUE})
#' to the model or not (\code{intercept = FALSE}).
#' 
#' @param tol The default tolerance is \code{epstol = 1e-12} and \code{convtol = 1e-8}.
#' See the documentation for \code{set_default_tolerance}. \code{epstol} stands for the
#' safeguard tolerance for mixture proportions (e.g. when \code{pi[1] * log(pi[1])} is
#' computed), and \code{convtol} stands for convergence tolerance.
#' 
#' @return A list object with the following elements:
#' 
#' \item{intercept}{An intercept.}
#' 
#' \item{beta}{A vector of estimated regression coefficients (variational posterior means), 
#' after fixed effects (e.g. intercept) from covariates are subtracted out.}
#' 
#' \item{sigma2}{A scalar value of estimated noise variance (approximate maximum likelihood).}
#' 
#' \item{pi}{A vector of estimated mixture proportions of length K, where \code{K = length(sa2)}.}
#' 
#' \item{iter}{The number of total outer loop iterations implemented in the coordinate ascent algorithm.}
#' 
#' \item{varobj}{A sequence of variational objective values (which equals the negative evidence lower bound).
#' \code{length(varobj)} should be equal to \code{iter}.}
#' 
#' \item{data}{A preprocessed data used as the actual input for the algorithm. When \code{Z = NULL}
#' and \code{intercept = TRUE}, then the columns of X and y will be centered, and returned.
#' In general, Z will be regressed out, or equivalently, X and y will be projected into the space
#' orthogonal to Z, and then will be returned.}
#' 
#' \item{update.order}{An update order used for the outer loop iterations.}
#' 
#' @references
#' 
#' Y. Kim, W. Wang, P. Carbonetto, M. Stephens (2020), Fast and
#'     Flexible Empirical Bayes Approach to Prediction in Multiple
#'     Regression.
#' 
#' @useDynLib mr.ash.alpha
#' 
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' @importFrom stats var
#' 
#' @examples
#' ### generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' ### fit Mr.ASH
#' fit.mr.ash  = mr.ash(X,y, method = "caisa")
#' 
#' ### prediction routine
#' Xnew        = matrix(rnorm(n*p),n,p)
#' ynew        = Xnew %*% beta + rnorm(n)
#' ypred       = predict(fit.mr.ash, Xnew)
#'
#' ### test error
#' rmse        = norm(ynew - ypred, '2') / sqrt(n)
#' 
#' ### coefficients
#' betahat     = predict(fit.mr.ash, type = "coefficients")
#' # this equals c(fit.mr.ash$intercept, fit.mr.ash$beta)
#' 
#' @export
#' 
mr.ash                      = function(X, y, Z = NULL, sa2 = NULL,
                                       method = c("caisa","sigma","accelerate",
                                                  "block","sigma_scaled",
                                                  "sigma_indep"),
                                       max.iter = 1000, min.iter = 1,
                                       beta.init = NULL,
                                       update.pi = TRUE, pi = NULL,
                                       update.sigma2 = TRUE, sigma2 = NULL,
                                       update.order = NULL,
                                       standardize = FALSE, intercept = TRUE,
                                       tol = set_default_tolerance()){
  
  # get sizes
  n            = nrow(X)
  p            = ncol(X)
  
  # check necessary conditions
  if (is.null(sa2)) {
    
  } else {
    if (any(sa2 < 0)) {
      stop ("all the mixture component variances must be non-negative.")
    }
    if (sa2[1] != 0) {
      stop ("the first mixture component variance sa2[1] must be 0.")
    }
  }
  
  # match method
  method       = match.arg(method)
  
  # set default tolerances unless specified
  tol0         = set_default_tolerance()
  tol          = modifyList(tol0,tol,keep.null = TRUE)
  
  # remove covariates
  data         = remove_covariate(X, y, Z, standardize, intercept)
  
  # initialize beta
  if ( is.null(beta.init) ){
    data$beta  = as.vector(double(p))
  } else {
    if (standardize) {
      data$beta  = as.vector(beta.init) * attr(data$X,"scaled:scale")
    } else {
      data$beta  = as.vector(beta.init)
    }
  }
  data$beta[1] = data$beta[1] + 0   # to make sure beta.init is not modified
  
  # initialize r
  r            = data$y - data$X %*% data$beta
  
  # sigma2
  if (is.null(sigma2))
    sigma2 = c(var(r))
  
  # set sa2 if missing
  if ( is.null(sa2) ) {
    if ( is.null(beta.init) )
      sa2      = (2^((0:19) / 20) - 1)^2
    else
      sa2      = (2^((0:19) / 5) - 1)^2
  }
  K            = length(sa2)
  data$sa2     = sa2
  
  # precalculate
  w            = colSums(data$X^2)
  data$w       = w
  
  # initialize other parameters
  if ( is.null(pi) ) {
    if ( is.null(beta.init) ){
      
      Phi          = matrix(1,p,K)/K
      pi           = rep(1,K)/K
      
    } else {
      
      S            = outer(1/w, sa2, '+') * sigma2
      Phi          = -data$beta^2/S/2 - log(S)/2
      Phi          = exp(Phi - apply(Phi,1,max))
      Phi          = Phi / rowSums(Phi)
      pi           = colMeans(Phi)
      
    }
  } else
    Phi          = matrix(rep(pi, each = p), nrow = p)
  
  # verbose = TRUE (TO DO LIST)
  verbose        = TRUE
  
  # run algorithm
  
  if ( is.null(update.order) ) {
    update.order   = 1:p
    if (method == "caisa") {
      if (update.pi) {
        out          = caisa_em   (data$X, w, sa2, pi, data$beta, r, sigma2,
                                   max.iter, min.iter, tol$convtol, tol$epstol,
                                   update.sigma2, verbose)
      } else {
        out          = caisa_fix_pi(data$X, w, sa2, pi, data$beta, r, sigma2,
                                    max.iter, min.iter, tol$convtol, tol$epstol,
                                    update.sigma2, verbose)
      }
    } else if (method == "sigma") {
      out          = caisa_sigma2  (data$X, w, sa2, pi, data$beta, r, sigma2,
                                    max.iter, min.iter, tol$convtol, tol$epstol,
                                    update.sigma2, verbose)
    } else if (method == "accelerate") {
      out           = caisa_acc (data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, mixsqpiter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
    } else if (method == "block") {
      stepsize      = 1
      out           = caisa_g  (data$X, w, sa2, Phi, pi, data$beta, r, sigma2,
                                max.iter, min.iter, tol$convtol, tol$epstol,
                                stepsize, update.sigma2, mode, verbose)
    } else if (method == "sigma_scaled") {
      out          = caisa_em2  (data$y, data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
      out$beta     = out$beta * sqrt(out$sigma2)
    } else if (method == "sigma_indep") {
      out          = caisa_em3  (data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
    }
  } else if (is.numeric(update.order)) {
    o   = rep(update.order - 1, max.iter)
    out = caisa_order(data$X, w, sa2, pi, data$beta, r, sigma2, 
                      o, max.iter, min.iter, tol$convtol, tol$epstol, 
                      update.sigma2, verbose)
  } else if (update.order == "random") {
    update.order = random_order(p, max.iter)
    out = caisa_order(data$X, w, sa2, pi, data$beta, r, sigma2, 
                      update.order, max.iter, min.iter, tol$convtol, tol$epstol, 
                      update.sigma2, verbose)
  }
  out$intercept     = c(data$ZtZiZy - data$ZtZiZX %*% out$beta)
  out$data         = data
  out$update.order = update.order
  if (standardize)
    out$beta = out$beta/attr(data$X, "scaled:scale")
  class(out) <- c("mr.ash", "list")
  
  if (out$pi[K] > tol$epstol) {
    warning("non-zero mixture proportion pi[K] on component with a largest variance sa2[K]: consider increasing sa2")
  }
  
  return(out)
}
