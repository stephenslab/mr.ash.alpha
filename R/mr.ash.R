#' @title Multiple Regression with Adaptive Shrinkage
#' 
#' @description Model fitting algorithms for Multiple Regression with
#'   Adaptive Shrinkage ("Mr.ASH"). Mr.ASH is a variational empirical
#'   Bayes (VEB) method for multiple linear regression. The fitting
#'   algorithms maximize the approximate marginal likelihood (the "evidence
#'   lower bound", or "ELBO") via coordinate-wise updates.
#' 
#' @details Mr.ASH adopts the following multiple linear regression
#'   model: \deqn{y | X, \beta, \sigma^2 \sim N(X \beta, \sigma^2 I_n),}
#'   in which the regression coefficients admit the following
#'   mixture-of-normals prior: \deqn{\beta | \pi, \sigma^2 ~
#'   \sum_{k=1}^K N(0, \sigma^2 \ sigma_k^2)}. Each mixture component is
#'   a normal density with zero mean and variance \eqn{\sigma^2
#'   \sigma_k^2}.
#' 
#'   The VEB approach solves the following optimization problem:
#'   \deqn{F(q,g,\sigma^2) = E_q \log p(y|X,\beta,\sigma^2) -
#'   \sum_{j=1}^p D_{KL}(q_j || g)} The algorithm updates the
#'   variational factors \eqn{q_1,...,q_p}, \eqn{g} and \eqn{\sigma^2}
#'   one at a time while fixing the others, in each outer loop
#'   iteration.
#' 
#'   See \sQuote{References} for more details about the VEB approach.
#' 
#' @seealso \code{\link{get.full.posterior}}
#'
#' @param X The input matrix, of dimension (n,p); each column is a
#'   single predictor; and each row is an observation vector. Here, n is
#'   the number of samples and p is the number of predictors. The matrix
#'   cannot be sparse.
#' 
#' @param y The observed quantitative responses, a vector of length p.
#' 
#' @param Z The covariate matrix, of dimension (n,k), where k is the
#'   number of covariates.  The input matrix Z can be modified according
#'   to "intercept" argument. If \code{Z = NULL} and \code{intercept =
#'   TRUE}, then the actual \code{Z} will be the matrix having entries 1
#'   of dimension (n,1). If \code{Z = NULL} and \code{intercept =
#'   FALSE}, no intercept or covariates are inclued the model. If
#'   \code{Z} is not \code{NULL} and \code{intercept = TRUE}, then the
#'   intercept is added as a covariate to \code{Z}.
#' 
#' @param sa2 The vector of mixture component variances. The variances
#'   should be in increasing order, starting at zero; that is,
#'   \code{sort(sa2)} should be the same as \code{sa2}. When \code{sa2 =
#'   NULL}, the default setting is used, \code{sa2[k] = (2^(0.05*(k-1))
#'   - 1)^2}, for \code{k = 1:20}. For this default setting,
#'   \code{sa2[1] = 0}, and \code{sa2[20]} is roughly 1.
#' 
#' @param method In the manuscript (see \sQuote{References}), only
#' \code{method = "caisa"} is used ("Cooridinate Ascent Iterative
#' Shinkage Algorithm"). Other method arguments will work, and produce
#' similar outcomes unless the regression setting is extreme.
#' (For dev 1) All the other settings but \code{method = "caisa"} are
#' purely experimental and under construction. 
#' The \code{method} arguments "block" and "accelerate" use different
#' updates for \eqn{g}.
#' (For dev 2) The \code{method} arguments "caisa", "sigma", "sigma_scaled",
#' "sigma_indep" use different updates for \eqn{sigma^2}, based on
#' different parametrizations on the variational posterior \eqn{q} and
#' \eqn{g}. More precisely, the update for \eqn{\sigma^2} depends on
#' whether we use sigma-dependent parametrization for \eqn{q} and/or
#' \eqn{g}. See reference for details. 
#' 
#' @param max.iter The maximum number of outer loop iterations allowed.
#' 
#' @param min.iter The minimum number of outer loop iterations allowed.
#' 
#' @param beta.init The initial estimate of the (approximate)
#'   posterior mean regression coefficients. This should be \code{NULL},
#'   or a vector of length p. When \code{beta.init = NULL}, the
#'   posterior mean coefficients are all initially zero.
#' 
#' @param update.pi If \code{update.pi = TRUE}, the mixture
#'   proportions in the mixture-of-normals prior are estimated from the
#'   data. In the manuscript, \code{update.pi = TRUE}.
#' 
#' @param pi The initial estimate of the mixture proportions
#'   \eqn{\pi_1,...,\pi_K}. If \code{pi = NULL}, the default value
#'   \code{pi[k] = 1/K} for \code{k = 1:K} will be used. 
#' 
#' @param update.sigma2 If \code{update.sigma2 = TRUE}, the residual
#'   variance \eqn{sigma^2} is estimated from the data.  In the manuscript,
#'   \code{update.sigma = TRUE}.
#' 
#' @param sigma2 The initial estimate of the residual variance,
#'   \eqn{\sigma^2}. If \code{sigma2 = NULL}, the residual variance is
#'   initialized to the empirical variance of the residuals based on the
#'   initial estimates of the regression coefficients, \code{beta.init},
#'   after removing linear effects of the intercept and any covariances.
#'
#' @param update.order The order in which the co-ordinate ascent
#'   updates for estimating the posterior mean coefficients are
#'   performed. \code{update.order} can be \code{NULL}, \code{"random"},
#'   or any permutation of \eqn{(1,...,p)}, where \code{p} is the number
#'   of columns in the input matrix \code{X}. When \code{update.order}
#'   is \code{NULL}, the co-ordinate ascent updates are performed in
#'   order in which they appear in \code{X}; this is equivalent to
#'   setting \code{update.order = 1:p}. When \code{update.order =
#'   "random"}, the co-ordinate ascent updates are performed in a
#'   randomly generated order.
#' 
#' @param standardize The logical flag for standardization of the
#'   columns of X variable, prior to the model fitting. The coefficients
#'   are always returned on the original scale.
#' 
#' @param intercept When \code{intercept = TRUE}, an intercept is
#'   included in the regression model.
#' 
#' @param tol Additional settings controlling behaviour of the model
#'   fitting algorithm. \code{tol$convtol} controls the termination
#'   criterion for the model fitting. When \code{update.pi = TRUE}, the
#'   outer-loop updates stop when the largest change in the mixture
#'   weights is less than \code{convtol*K}; when \code{update.pi =
#'   FALSE}, the outer-loop updates stop when the largest change in the
#'   estimates of the posterior mean coefficients is less than
#'   \code{convtol*K}. \code{tol$epstol} is a small, positive number
#'   added to the likelihoods to avoid logarithms of zero.
#' 
#' @return A list object with the following elements:
#' 
#' \item{intercept}{The estimated intercept.}
#' 
#' \item{beta}{A vector containing posterior mean estimates of the
#'   regression coefficients for all predictors.}
#' 
#' \item{sigma2}{The estimated residual variance.}
#' 
#' \item{pi}{A vector of containing the estimated mixture
#'   proportions}.
#' 
#' \item{iter}{The number of outer-loop iterations that were
#'   performed.}
#' 
#' \item{update.order}{The ordering used for performing the
#'   coordinate-wise updates.}
#' 
#' \item{varobj}{A vector, with \code{length(varobj) = numiter},
#'   containing the value of the variational objective (equal to the
#'   negative of the "evidence lower bound") attained at each outer-loop
#'   iteration of the model fitting algorithm.}
#' 
#' \item{data}{A preprocessed data used as the actual input for the
#' algorithm. When \code{Z = NULL} and \code{intercept = TRUE}, then
#' the columns of X and y will be centered, and returned.  In general,
#' Z will be regressed out, or equivalently, X and y will be projected
#' into the space orthogonal to Z, and then will be returned.}
#' 
#' @seealso \code{\link{get.full.posterior}}
#' 
#' @references
#'
#' Y. Kim, W. Wang, P. Carbonetto, M. Stephens (2020). Fast and
#' flexible empirical Bayes approach to prediction in multiple
#' regression.
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
                                       method = c("caisa","sigma","sigma_scaled",
                                                  "sigma_indep","accelerate",
                                                  "block"),
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
  if (!is.null(sa2)) {
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
    sa2      = (2^((0:19) / 20) - 1)^2
  }
  K            = length(sa2)
  data$sa2     = sa2
  
  # precompute x_j^T x_j
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
    } else if (method == "sigma_scaled") {
      out          = caisa_em2  (data$y, data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
      out$beta     = out$beta * sqrt(out$sigma2)
    } else if (method == "sigma_indep") {
      out          = caisa_em3  (data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
    } else if (method == "accelerate") {
      mixsqpiter    = 5
      out           = caisa_acc (data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, mixsqpiter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
    } else if (method == "block") {
      stepsize      = 1
      out           = caisa_g  (data$X, w, sa2, Phi, pi, data$beta, r, sigma2,
                                max.iter, min.iter, tol$convtol, tol$epstol,
                                stepsize, update.sigma2, mode, verbose)
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
    warning("The mixture proportion associated with the largest prior variance ",
            "is greater than zero; this indicates that the model fit could be ",
            "improved by using a larger setting of the prior variance. Consider ",
            "increasing the range of the variances \"sa2\".")
  }
  
  return(out)
}
