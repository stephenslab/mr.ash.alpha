# ----------------------------------------------------------------------
# Some miscellaneuous auxiliary functions are listed above.
# Some functions are directly copied from varbvs,
# https://github.com/pcarbo/varbvs
# ----------------------------------------------------------------------

# Remove covariate effects Regresses Z out from X and y; that is, X
# and y are projected into the space orthogonal to Z.
#' 
#' @importFrom Matrix forceSymmetric
#'
remove_covariate <- function (X, y, Z, standardize = FALSE, intercept = TRUE) {
  
  # check if Z is null and intercept = FALSE
  if (is.null(Z) & (intercept == FALSE)) {
    return(list(X = X, y = y, Z = Z,
                ZtZiZX = rep(0,dim(X)[2]), ZtZiZy = 0))
  }
  
  # redefine y
  y = c(as.double(y))
  n = length(y)
  
  # add intercept if intercept = TRUE
  if (intercept) {
    if (is.null(Z))
      Z <- matrix(1,n,1)
    else
      Z <- cbind(1,Z)
  }
  
  if (ncol(Z) == 1) {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    X           = scale(X, center = intercept, scale = standardize)
    alpha       = mean(y)
    y           = y - alpha
    
  } else {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    
    #   y = y - Z (Z^T Z)^{-1} Z^T y
    #   X = X - Z (Z^T Z)^{-1} Z^T X  
    y     = y - c(Z %*% ZtZiZy)
    X     = X - Z %*% ZtZiZX
  }
  
  return(list(X = X, y = y, Z = Z,
              ZtZiZX = ZtZiZX, ZtZiZy = ZtZiZy))
}

#' @title Ordering of Predictors from Univariate Regression
#' 
#' @description This function extracts the ordering of the predictors
#'   according to the coefficients estimated in a basic univariate
#'   regression; in particular, the predictors are ordered in decreasing
#'   order by magnitude of the univariate regression coefficient
#'   estimate.
#' 
#' @param X An input design matrix. This may be centered and/or
#'   standardized prior to calling function.
#' 
#' @param y A vector of response variables.
#'
#' @return An ordering of the predictors.
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
#' univ.order = univar.order(X,y)
#' 
#' @export
#' 
univar.order = function(X, y) {
  colnorm = c(colMeans(X^2))
  return (order(abs(c(t(X) %*% y) / colnorm), decreasing = TRUE))
}

#' @title Ordering of Predictors from Coefficient Estimates 
#' 
#' @param beta A vector of estimated regression coefficients.
#' 
#' @description This function orders the predictors by decreasing
#'   order of the magnitude of the estimated regression coefficient.
#'
#' @return An ordering of the predictors.
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
#' ### glmnet fit
#' library(glmnet)
#' beta.lasso = coef(cv.glmnet(X, y))[-1]
#' lasso.order = absolute.order(beta.lasso)
#' 
#' ### ncvreg fit
#' library(ncvreg)
#' beta.scad = c(coef(cv.ncvreg(X, y))[-1])
#' scad.order = absolute.order(beta.scad)
#' 
#' @export
#' 
absolute.order = function (beta) {
  abs_order = c(order(abs(beta), decreasing = TRUE))
  return (abs_order)
}

#' @title Ordering of Predictors by Regularization Path
#' 
#' @param fit The output of a function such as \code{glmnet} from the
#'   \code{glmnet} package or \code{ncvreg} from the \code{ncvfeg} that
#'   estimates a "regularization path" for all predictors.
#' 
#' @description This function determines an ordering of the predictors
#'  based on the regularization path of the penalized regression; in
#'   particular, the predictors are ordered based on the order in which
#'   the coefficients are included in the model as the penalty strength
#'   decreases.
#' 
#' @return An ordering of the predictors.
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
#' ### glmnet fit
#' library(glmnet)
#' fit.lasso = glmnet(X, y)
#' lasso.order = path.order(fit.lasso)
#' 
#' ### ncvreg fit
#' library(ncvreg)
#' fit.scad = ncvreg(X, y)
#' scad.order = path.order(fit.scad)
#'
#' @export
#' 
path.order = function (fit) {
  beta_path = coef(fit)[-1,]
  K = dim(beta_path)[2]
  path_order = c()
  for (k in 1:K) {
    crt_path = which(beta_path[,k] != 0)
    if (length(crt_path) != 0 & length(path_order) == 0) {
      path_order = c(path_order, crt_path)
    } else if(length(crt_path) != 0) {
      path_order = c(path_order, crt_path[-which(crt_path %in% path_order)] )
    }
  }
  path_order = unname(path_order)
  index_order = c(path_order, seq(1,dim(beta_path)[1])[-path_order])
  return (index_order)
}

#' @title Extract Regression Coefficients from Mr.ASH Fit
#'
#' @description Retrieve posterior mean estimates of the regression
#'   coefficients in a Mr.ASH model.
#' 
#' @param object A Mr.ASH fit, usually the result of calling
#'   \code{mr.ash}.
#'
#' @param ... Additional arguments passed to the default S3 method.
#' 
#' @return A p+1 vector. The first element gives the estimated
#'   intercept, and the remaining p elements are the estimated
#'   regression coefficients.
#'   
#' ## generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' ## fit mr.ash
#' fit.mr.ash  = mr.ash(X, y)
#' 
#' ## coefficient
#' coef.mr.ash = coef(fit.mr.ash)
#' intercept   = coef.mr.ash[1]
#' beta        = coef.mr.ash[-1]
#' 
#' @importFrom stats coef
#' 
#' @export coef.mr.ash
#' 
#' @export
#' 
coef.mr.ash = function (object, ...)
  c(object$intercept,object$beta)

#' @title Predict Outcomes or Extract Coefficients from Mr.ASH Fit
#'
#' @description This function predicts outcomes (y) given the observed
#'   variables (X) and a Mr.ASH model; alternatively, retrieve the
#'   estimates of the regression coefficients.
#'
#' @param object A mr_ash fit, usually the result of calling
#'   \code{mr.ash}.
#'
#' @param newx The input matrix, of dimension (n,p); each column is a
#'   single predictor; and each row is an observation vector. Here, n is
#'   the number of samples and p is the number of predictors. When
#'   \code{newx} is \code{NULL}, the fitted values for the training data
#'   are provided.
#' 
#' @param type The type of output. For \code{type = "response"},
#'   predicted or fitted outcomes are returned; for \code{type =
#'   "coefficients"}, the estimated coefficients are returned.
#' 
#' @param ... Additional arguments passed to the default S3 method.
#'
#' @return For \code{type = "response"}, predicted or fitted outcomes
#' are returned; for \code{type = "coefficients"}, the estimated
#' coefficients are returned.
#' 
#' @examples
#' ## generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' ## fit mr.ash
#' fit.mr.ash  = mr.ash(X, y)
#' 
#' ## predict
#' Xnew        = matrix(rnorm(n*p),n,p)
#' ypred       = predict(fit.mr.ash, Xnew)
#' 
#' @importFrom stats predict
#' 
#' @export predict.mr.ash
#' 
#' @export
#' 
predict.mr.ash               = function(object,newx = NULL,
                                        type=c("response","coefficients"),...) {
  
  type <- match.arg(type)
  if (type == "coefficients"){
    if(!missing(newx))
      stop("Do not supply newx when predicting coefficients")
    return(coef(object))
  }
  else if(missing(newx))
    return(object$fitted)
  else {
    if (!all(object$data$Z == 1))
      stop("predict.mr.ash is not implemented for covariates Z other than ",
           "intercept")
    return(drop(object$intercept + newx %*% coef(object)[-1]))
  }
}

set_default_tolerance       = function(){
  epstol    = 1e-12
  convtol   = 1e-8
  
  return ( list(epstol = epstol, convtol = convtol ) )
}

#' @title Approximation Posterior Expectations from Mr.ASH Fit
#'
#' @description Recover the parameters specifying the variational
#'   approximation to the posterior distribution of the regression
#'   coefficients. To streamline the model fitting implementation, and
#'   to reduce memory requirements, \code{\link{mr.ash}} does not store
#'   all the parameters needed to specify the approximate posterior.
#' 
#' @param fit A Mr.ASH fit obtained, for example, by running
#'   \code{mr.ash}.
#' 
#' @return A list object with the following elements:
#' 
#' \item{phi}{An p x K matrix containing the posterior assignment
#'   probabilities, where p is the number of predictors, and K is the
#'   number of mixture components. (Each row of \code{phi} should sum to
#'   1.)}
#' 
#' \item{m}{An p x K matrix containing the posterior means conditional
#'   on assignment to each mixture component.}
#' 
#' \item{s2}{An p x K matrix containing the posterior variances
#'   conditional on assignment to each mixture component.}
#' 
#' @examples
#' ## generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' ## fit mr.ash
#' fit.mr.ash  = mr.ash(X, y)
#' 
#' ## recover full posterior
#' full.post   = get.full.posterior(fit.mr.ash)
#' 
#' @export
#' 
get.full.posterior <- function(fit) {
    
  # compute residual
  r = fit$data$y - fit$data$X %*% fit$beta
  
  # compute bw and s2
  bw = as.vector((t(fit$data$X) %*% r) + fit$data$w * fit$beta)
  s2 = fit$sigma2 / outer(fit$data$w, 1/fit$data$sa2, '+')
  
  # compute m, phi
  m   = bw * s2 / fit$sigma2
  phi = -log(1 + outer(fit$data$w,fit$data$sa2))/2 + m * (bw/2/fit$sigma2)
  phi = c(fit$pi) * t(exp(phi - apply(phi,1,max)))
  phi = t(phi) / colSums(phi)
  return (list(phi = phi, m = m, s2 = s2))
}

gibbs.sampling              = function(X, y, pi, sa2 = (2^((0:19) / 20) - 1)^2,
                                       max.iter = 1500, burn.in = 500,
                                       standardize = FALSE, intercept = TRUE,
                                       sigma2 = NULL, beta.init = NULL,
                                       verbose = TRUE){
  
  # get sizes
  n            = nrow(X)
  p            = ncol(X)
  
  # remove covariates
  data         = remove_covariate(X, y, NULL, standardize, intercept)
  if ( is.null(beta.init) )
    data$beta  = as.vector(double(p))
  else
    data$beta  = as.vector(beta.init)
  
  # initialize r
  r            = data$y - data$X %*% data$beta
  
  # sigma2
  if ( is.null(sigma2) )
    sigma2 = c(var(r))
  
  # precalculate
  w            = colSums(data$X^2)
  data$w       = w
  
  # gibbs sampling
  out           = gibbs.sampling(data$X, w, sa2, pi, data$beta, r, sigma2, max.iter, burn.in, verbose)
  out$data      = data
  out$mu        = c(data$ZtZiZy - data$ZtZiZX %*% out$beta)
  
  return (out)
}

var.n                       = function(x) {
  a             = x - mean(x)
  return (sum(a^2) / length(a))
}
