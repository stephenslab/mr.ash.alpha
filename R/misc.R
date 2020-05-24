# ----------------------------------------------------------------------
# Some miscellaneuous auxiliary functions are listed above.
# Some functions are directly borrowed from the R package varbvs
# https://github.com/pcarbo/varbvs
# ----------------------------------------------------------------------

#' @title remove covariate effects
#' @description This function regresses the effect of Z out from X and y.
#' In other words, X and y will be projected into the space orthogonal to Z.
#' 
#' @importFrom Matrix forceSymmetric
#' 
remove_covariate <- function (X, y, Z, standardize = FALSE, intercept = TRUE) {
  
  # redefine y
  
  y = c(as.double(y))
  n = length(y)
  
  # add intercept if intercept == TRUE
  
  if (intercept == TRUE) {
    if (is.null(Z))
      Z <- matrix(1,n,1)
    else
      Z <- cbind(1,Z)
  } else {
    return(list(X = X, y = y))
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
    alpha = c(Z %*% ZtZiZy)
    y     = y - alpha
    X     = X - Z %*% ZtZiZX
  }
  
  return(list(X = X, y = y, Z = Z, alpha = alpha, ZtZiZX = ZtZiZX, ZtZiZy = ZtZiZy))
}

#' @title regularization path order
#' 
#' @param beta a glmnet fit, or a ncvreg fit
#' 
#' @description This function extracts the path order from penalized regression fit.
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
#' lasso.order = mr.ash.alpha:::path.order(fit.lasso)
#' 
#' ### ncvreg fit
#' library(ncvreg)
#' fit.scad = ncvreg(X, y)
#' scad.order = mr.ash.alpha:::path.order(fit.scad)
#' 
#' @export
#' 
path.order = function (fit) {
  # perform lasso regression and reorder regressors by "importance"
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

#' @title absolute magnitude order
#' 
#' @param beta a vector of regression coefficients
#' 
#' @description This function extracts the absolute value order from any fit.
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
#' lasso.order = mr.ash.alpha:::abs.order(beta.lasso)
#' 
#' ### ncvreg fit
#' library(ncvreg)
#' beta.scad = c(coef(cv.ncvreg(X, y))[-1])
#' scad.order = mr.ash.alpha:::abs.order(beta.scad)
#' 
#' @export
#' 
absolute.order = function (beta) {
  
  # abs order
  abs_order = c(order(abs(beta), decreasing = TRUE))
  return (abs_order)
}

#' @title univariate order
#' 
#' @param X An input design matrix.
#' @param y A vector of response variables.
#' 
#' @description This function extracts the univariate regression coefficient order
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

#' @title extract regression coefficients from mr_ash fit
#' 
#' @param object a mr_ash fit
#' 
#' @return a p+1 vector, the first element being an intercept, and the
#'   remaining p elements being estimated regression coefficients
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
#' @export coef.mr.ash
#' @export
#' 
coef.mr.ash = function (object, ...)
  c(object$intercept,object$beta)

#' @title predict future observations or extract coefficients from mr_ash fit
#' 
#' @param object a mr_ash fit
#' 
#' @param newx a new value for X at which to do predictions
#' 
#' @param type if this is coefficients, then calls coef.susie
#' 
#' @importFrom stats coef
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
#' 
#' @export predict.mr.ash
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
  
  if(missing(newx)){return(object$fitted)}
  
  return(drop(object$intercept + newx %*% coef(object)[-1]))
}

#' @title set_default_tolerance
#' 
#' @description 
#' 
set_default_tolerance       = function(){
  epstol    = 1e-12
  convtol   = 1e-8 #sqrt(.Machine$double.eps)
  
  return ( list(epstol = epstol, convtol = convtol ) )
}

#' @title get q
#' 
#' @param fit An Mr.ASH fit obtained by running \code{mr.ash}. See \sQuote{Examples}.
#' 
#' @return 
#' A list object with the following elements:
#' 
#' \item{phi}{A posterior component probability matrix of dimension (p,K);
#' each row is a vector posterior mixture proportions corresponding to
#' each regression coefficient.}
#' 
#' \item{m}{A posterior component mean matrix of dimension (p,K);
#' each row is a vector of posterior component means corresponding to
#' each regression coefficient.}
#' 
#' \item{s2}{A posterior component variance matrix of dimension (p,K);
#' each row is a vector of posterior component variances corresponding to
#' each regression coefficient.}
#' 
#' @description The \code{get.full.posterior} function recovers the estimated
#' variational posterior, which maximizes the ELBO F(q,g,sigma2). The \code{mr.ash}
#' function does not store the full posterior for a faster implementation, thus
#' this function \code{get.full.posterior} recovers the full posterior from the
#' Mr.ASH fit.
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
  r            = fit$data$y - fit$data$X %*% fit$beta
  
  # compute bw and s2
  bw           = as.vector((t(fit$data$X) %*% r) + fit$data$w * fit$beta)
  s2           = 1 / outer(fit$data$w, 1/fit$data$sa2, '+');
  
  # compute m, phi
  m            = bw * s2;
  phi          = -log(1 + outer(fit$data$w, fit$data$sa2))/2 + m * (bw / 2 / fit$sigma2);
  phi          = c(fit$pi) * t(exp(phi - apply(phi,1,max)));
  phi          = t(phi) / colSums(phi);
  return (list(phi = phi, m = m, s2 = s2))
}


#' @title gibbs sampling
#' @export
#' 
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
  out           = gibbs_sampling(data$X, w, sa2, pi, data$beta, r, sigma2, max.iter, burn.in, verbose)
  out$data      = data
  out$mu        = c(data$ZtZiZy - data$ZtZiZX %*% out$beta)
  
  return (out)
}
