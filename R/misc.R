#' ----------------------------------------------------------------------
#' Some miscellaneuous auxiliary functions are listed above.
#' Some functions are directly borrowed from the R package varbvs
#' https://github.com/pcarbo/varbvs
#' ----------------------------------------------------------------------

#' ----------------------------------------------------------------------
#' y = Z * alpha + X * beta + epsilon
#' remove Z and project X and y onto space orthogonal to Z
#' add intercept (a vecror of ones) to Z if intercept == TRUE
#' standardize X if standardize == TRUE
#' ----------------------------------------------------------------------
#' @title remove covariate effects
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

#' @export
# define order function
path.order = function (fit.glmnet) {
  # perform lasso regression and reorder regressors by "importance"
  beta_path = coef(fit.glmnet)[-1,]
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

#' @export
abs.order = function(beta) {
  return (order(abs(beta), decreasing = TRUE))
}

#' @export
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

#' ----------------------------------------------------------------------
#'
#'
#' ----------------------------------------------------------------------
set_default_tolerance       = function(){
  epstol    = 1e-10
  convtol   = 1e-6 #sqrt(.Machine$double.eps)
  
  return ( list(epstol = epstol, convtol = convtol ) )
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
