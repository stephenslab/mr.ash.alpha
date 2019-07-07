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

#' @title extract regression coefficients from mr_ash fit
#' @param object a mr_ash fit
#' @return a p+1 vector, the first element being an intercept, and the remaining p elements being estimated regression coefficients
#' @export coef.mr_ash
#' @export
coef.mr_ash = function(object, ...){
  c(object$intercept, object$beta)
}

#' @title predict future observations or extract coefficients from mr_ash fit
#' @param object a mr_ash fit
#' @param newx a new value for X at which to do predictions
#' @param type if this is coefficients, then calls coef.susie
#' @importFrom stats coef
#' @export predict.mr_ash
#' @export
#' @export
predict.mr_ash               = function(object,newx = NULL,
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


sim_data = function(n, p, pi0 = 0.85, signal.shape = "spiky", center = TRUE,
                    corr = 0, seed = 2018) {
  
  # set seed
  set.seed(seed);
  
  # construct covariance structure of X
  Sigma = matrix(0,p,p);
  if (corr > 0) {
    for(i in 1:p){
      Sigma[i:p,i] = exp(-((0:(p-i)))^2/corr);
      Sigma[i,i:p] = exp(-((0:(p-i)))^2/corr);
    }
    L         = chol(Sigma)                        # Sigma = t(L) %*% L
    X         = matrix(rnorm(n*p),nrow = n) %*% L  # X has a column covariance Sigma in population level
  } else if (corr == 0) {
    X = matrix(rnorm(n*p),nrow = n);                 # X has a column covariance I_p in population level
  } else {
    stop("corr must be greater than 0");
  }
  
  # construct pi
  pi        = c(pi0, rep((1-pi0)/3, 3))
  pi        = pi / sum(pi)
  
  # construct label
  label     = sample(length(pi), p, replace = TRUE, prob = pi)
  
  
  # construct signal
  if (signal.shape == "spiky") {
    sd        = c(0,1,2,6);
    signal    = 0.45;
    beta      = signal * rnorm(p, 0, sd[label]);
  } else if (signal.shape == "nearnormal") {
    sd        = c(0,3,4,5);
    signal    = 0.4;
    beta      = signal * rnorm(p, 0, sd[label]);
  } else if (signal.shape == "flattop") {
    sd        = c(0,1.5,1.75,2);
    signal    = 0.5;
    beta      = c(signal * rnorm(p, 10 * (runif(p) - 1/2), sd[label]))
  } else if (signal.shape == "bimodal") {
    sd        = c(0,1,3,5);
    signal    = 0.4;
    beta      = c(signal * rnorm(p, 5 * rbinom(p,1,.5) - 2.5, sd[label]))
  } else if (signal.shape == "heavytail") {
    signal    = 6;
    beta      = c(signal * rexp(p, 5) * sign(rnorm(p)))
  } else if (signal.shape == "polygenic") {
    pi        = c(pi0, (1-pi0) * 0.495, (1-pi0) * 0.495, (1-pi0) * 0.01)
    label     = sample(length(pi), p, replace = TRUE, prob = pi)
    sd        = c(0,0.05,0.1,6);
    signal    = 0.5;
    beta      = signal * rnorm(p, 0, sd[label]);
  }
  beta[label == 1] = 0;
  
  if (center == TRUE) {
    X         = t(t(X) - colMeans(X))
    y         = X %*% beta + rnorm(n);
    y         = y - mean(y)
  } else {
    y         = X %*% beta + rnorm(n);
  }
  
  return ( list ( X = X, y = y, beta = beta, label = label, Sigma = Sigma[1:10,1:10]) )
}

generate_beta = function(p, pi0 = 0.95, signal.shape = "spiky", signal.strength = 1, seed = 2019) {
  
  # set.seed
  set.seed(seed)
  
  # construct pi
  pi        = c(pi0, rep((1-pi0)/3, 3))
  pi        = pi / sum(pi)
  
  # construct label
  label     = sample(length(pi), p, replace = TRUE, prob = pi)
  
  # construct signal
  if (signal.shape == "spiky") {
    sd        = c(0,1,2,6);
    signal    = 0.45;
    beta      = signal * rnorm(p, 0, sd[label]);
  } else if (signal.shape == "nearnormal") {
    sd        = c(0,3,4,5);
    signal    = 0.4;
    beta      = signal * rnorm(p, 0, sd[label]);
  } else if (signal.shape == "flattop") {
    sd        = c(0,1.5,1.75,2);
    signal    = 0.5;
    beta      = c(signal * rnorm(p, 10 * (runif(p) - 1/2), sd[label]))
  } else if (signal.shape == "bimodal") {
    sd        = c(0,1,3,5);
    signal    = 0.4;
    beta      = c(signal * rnorm(p, 5 * rbinom(p,1,.5) - 2.5, sd[label]))
  } else if (signal.shape == "heavytail") {
    signal    = 6;
    beta      = c(signal * rexp(p, 5) * sign(rnorm(p)))
  } else if (signal.shape == "polygenic") {
    pi        = c(pi0, (1-pi0) * 0.495, (1-pi0) * 0.495, (1-pi0) * 0.01)
    label     = sample(length(pi), p, replace = TRUE, prob = pi)
    sd        = c(0,0.05,0.1,6);
    signal    = 0.5;
    beta      = signal * rnorm(p, 0, sd[label]);
  }
  beta[label == 1] = 0;
  
  return (c(beta * signal.strength))
}


## compute posterior quantile

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