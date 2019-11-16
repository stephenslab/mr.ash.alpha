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

# define order function
path.order = function(fit.glmnet) {
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

abs.order = function(beta) {
  return (order(abs(beta), decreasing = TRUE))
}

univar.order = function(X, y) {
  colnorm = c(colMeans(X^2))
  return (order(abs(c(t(X) %*% y) / colnorm), decreasing = TRUE))
}

#' @title extract regression coefficients from mr_ash fit
#' @param object a mr_ash fit
#' @return a p+1 vector, the first element being an intercept, and the remaining p elements being estimated regression coefficients
#' @export coef.mr.ash
#' @export
coef.mr.ash = function(object, ...){
  c(object$intercept, object$beta)
}

#' @title predict future observations or extract coefficients from mr_ash fit
#' @param object a mr_ash fit
#' @param newx a new value for X at which to do predictions
#' @param type if this is coefficients, then calls coef.susie
#' @importFrom stats coef
#' @export predict.mr.ash
#' @export
predict.mr.ash               = function(object,newx = NULL,
                                        type=c("response","coefficients"),...) {
  
  type <- match.arg(type)
  if (type=="coefficients"){
    if(!missing(newx)){
      stop("Do not supply newx when predicting coefficients")
    }
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
  epstol    = 1e-10;
  convtol   = 1e-6; #sqrt(.Machine$double.eps);
  
  return ( list(epstol = epstol, convtol = convtol ) )
}


#' ----------------------------------------------------------------------
#'
#'
#' ----------------------------------------------------------------------
logspace                    = function(a = 0, b, length = 10){
  
  if (a >= b){
    stop("b must be strictly greater than a")
  }
  
  if (a > 0){
    return( exp( seq(log(a), log(b), length = length) ) )
  } else if (a == 0){
    return( exp( seq(0, log(b+1), length = length ) ) - 1 )
  } else{
    stop("a must be greater than 0")
  }
}


#' ----------------------------------------------------------------------
#' Compute the local false sign rate (LFSR) for each variable. This
#' assumes that the first mixture component is a "spike" (that is, a
#' normal density with a variance approaching zero).
#' ----------------------------------------------------------------------
#' @title compute local false discovery rate
#' @importFrom stats pnorm rnorm runif rbinom rexp
#' 
computelfsrmix <- function (alpha, mu, s) {
  
  # Get the number of variables (p) and the number of mixture
  # components (k).
  p <- nrow(alpha)
  k <- ncol(alpha)
  
  # For each variable, get the posterior probability that the
  # regression coefficient is exactly zero.
  p0 <- alpha[,1]
  
  # For each variable, get the posterior probability that the
  # regression coefficient is negative.
  if (k == 2)
    pn <- alpha[,2] * pnorm(0,mu[,2],sqrt(s[,2]))
  else
    pn <- rowSums(alpha[,-1] * pnorm(0,mu[,-1],sqrt(s[,-1])))
  
  # Compute the local false sign rate (LFSR) following the formula
  # given in the Biostatistics paper, "False discovery rates: a new
  # deal".
  lfsr     <- rep(0,p)
  b        <- pn > 0.5*(1 - p0)
  lfsr[b]  <- 1 - pn[b]
  lfsr[!b] <- p0[!b] + pn[!b]
  
  return(lfsr)
}

get_phi <- function(fit) {
  # compute residual
  r            = fit$data$y - fit$data$X %*% fit$beta
  
  # compute bw and S2inv
  bw           = as.vector((t(fit$data$X) %*% r) + fit$data$w * fit$beta)
  S2inv        = 1 / outer(fit$data$w, 1/fit$data$sa2, '+');
  
  # compute mu, phi
  mu           = bw * S2inv;
  phi          = -log(1 + outer(fit$data$w, fit$data$sa2))/2 + mu * (bw / 2 / fit$sigma2);
  phi          = c(fit$pi) * t(exp(phi - apply(phi,1,max)));
  phi          = t(phi) / colSums(phi);
  return (list(phi = phi, mu = mu, r = r))
}

#' @title compute posterior quantile
#' @importFrom stats pnorm uniroot
#' @export
#' 
posterior_quantile  <- function(fit, thresh = NULL) {
  
  # compute residual
  r            = fit$data$y - fit$data$X %*% fit$beta
  
  # compute bw and S2inv
  bw           = as.vector((t(fit$data$X) %*% r) + fit$data$w * fit$beta)
  S2inv        = 1 / outer(fit$data$w, 1/fit$data$sa2, '+');
  
  # compute mu, phi
  mu           = bw * S2inv;
  phi          = -log(1 + outer(fit$data$w, fit$data$sa2))/2 + mu * (bw / 2 / fit$sigma2);
  phi          = c(fit$pi) * t(exp(phi - apply(phi,1,max)));
  phi          = t(phi) / colSums(phi);
  
  if (is.null(thresh))
    return (phi)
  
  # calculate component posterior sds
  cpm          = mu[,-1]
  cpsd         = 1 / sqrt(outer(fit$data$w,1/fit$data$sa2[-1],'+'));
  
  # nonzero indices
  postquant    = matrix(0,length(bw),length(thresh));
  
  # compute posterior quantiles
  for (t in 1:length(thresh)) {
    
    ind          = which(phi[,1] < thresh[t]);
    
    for (i in ind) {
      
      # function for posterior median calculation
      # a unique zero of the function post_cdf is the posterior median
      post_cdf <- function(x){
        out = sum(phi[i,-1] * pnorm(x, abs(cpm[i,]), cpsd[i,])) + phi[i,1] - thresh[t]
      }
      
      if (post_cdf(0) >= 0) {
        postquant[i,t] = 0;
      } else {
        postquant[i,t] = uniroot(post_cdf, interval = c(0,1e8))$root * sign(bw[i]);
      }
    }
    
  }
  
  return (postquant)
}