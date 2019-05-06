#' ----------------------------------------------------------------------
#' Some miscellaneuous auxiliary functions are listed above.
#' Some functions are directly borrowed from the R package varbvs
#' https://github.com/pcarbo/varbvs
#' ----------------------------------------------------------------------


post_median      <- function(beta, r, d, Phi, sa2, thresh = 0.9) {
  
  # take step
  b            = as.vector((t(X) %*% r) / d + beta)
  
  # calculate component posterior means
  cpm          = b / (1 + 1 / outer(d,sa2[-1],'*'));
  
  # calculate component posterior sds
  cpsd         = 1 / sqrt(outer(d,1/sa2[-1],'+'));
  
  # nonzero indices
  ind          = which(Phi[,1] < thresh);
  postmed      = double(length(b));
  
  for (i in ind) {

    # function for posterior median calculation
    # a unique zero of the function post_cdf is the posterior median
    post_cdf <- function(x){
      out = sum(Phi[i,-1] * pnorm(x, abs(cpm[i,]), cpsd[i,])) + Phi[i,1] - thresh
    }
    
    if (post_cdf(0) >= 0) {
      postmed[i] = 0;
    } else {
      postmed[i] = uniroot(post_cdf, interval = c(0,1e3))$root * sign(b[i]);
    }
  }
  
  return (postmed)
}

#' ----------------------------------------------------------------------
#' y = Z * alpha + X * beta + epsilon
#' remove Z and project X and y onto space orthogonal to Z
#' add intercept (a vecror of ones) to Z if intercept == TRUE
#' standardize X if standardize == TRUE
#' ----------------------------------------------------------------------
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

#' ----------------------------------------------------------------------
#'
#'
#' ----------------------------------------------------------------------
predict.mr_ash               = function(fit, X){
  #scale(X, center = TRUE, scale = FALSE) %*% fit$beta + fit$data$alpha
  X %*% fit$beta + c(fit$intercept)
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
#'
#'
#' ----------------------------------------------------------------------
proj_l0                     = function(b, K = 10){
  
  # projection onto l0 ball
  b[order(b, decreasing = TRUE)[1:(length(b) - K)]] = 0
  
  return(b)
  
}


#' ----------------------------------------------------------------------
#' Compute the local false sign rate (LFSR) for each variable. This
#' assumes that the first mixture component is a "spike" (that is, a
#' normal density with a variance approaching zero).
#' ----------------------------------------------------------------------
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