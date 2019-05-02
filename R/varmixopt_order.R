varmixopt_order             = function(X, y, Z = NULL, sa2 = NULL,
                                       order.method = "1:p", o = NULL, K = 10,
                                       max.iter = 1000, min.iter = 1,
                                       stepsize = 1, 
                                       tol = list(), outputlevel = 1,
                                       update.sigma = TRUE, sigma2 = NULL,
                                       beta.init = NULL,
                                       standardize = FALSE, intercept = TRUE,
                                       verbose = TRUE){
  
  # stop if standardize = TRUE
  if ( standardize )
    stop("Currently we do not support standardization.")
  
  # get sizes
  n            = dim(X)[1];
  p            = dim(X)[2];
  
  # set default tolerances unless specified
  tol0         = set_default_tolerance()
  tol          = modifyList(tol0,tol,keep.null = TRUE)
  
  # calculate sa2 if missing
  if ( is.null(sa2) ) {
    betahat    = (t(X) %*% y) / n
    sa2        = logspace( b = 2 * ceiling(max(betahat^2)), length = K )
  }
  K            = length(sa2);
  
  # remove covariates
  data         = remove_covariate(X, y, Z, standardize, intercept);
  data$sa2     = sa2
  
  # initialize beta
  if ( is.null(beta.init) ){
    data$beta  = as.vector(double(p));
  } else {
    data$beta  = as.vector(beta.init);
  }
  data$beta[1] = data$beta[1] + 0;      # to make sure beta.init is not modified
  
  # initialize r
  r            = data$y - data$X %*% data$beta;
  
  # sigma2
  if ( is.null(sigma2) ) {
    sigma2 = c(var(r));
  }
  
  # precalculate
  w            = colSums(data$X^2);
  data$w       = w;
  
  # initialize other parameters
  if ( is.null(beta.init) ){
    #S            = outer(1/w, sa2, '+');
    Phi          = matrix(1,p,K)/K;
    pi           = rep(1,K)/K;
    #Phi          = -data$beta^2/S/2 - log(S)/2;
    #Phi          = exp(Phi);
    #Phi          = Phi / rowSums(Phi);
    #pi           = colMeans(Phi);
  } else {
    S            = outer(1/w, sa2, '+') * sigma2;
    Phi          = -data$beta^2/S/2 - log(S)/2;
    Phi          = exp(Phi - apply(Phi,1,max));
    Phi          = Phi / rowSums(Phi);
    pi           = colMeans(Phi);
  }
  
  if (order.method == "normal") {
    o = rep(0:(p-1), max.iter)
  } else if (order.method == "random") {
    o = sample(0:(p-1), p)
    for (i in 2:max.iter) {
      o = c(o, sample(0:(p-1), p))
    }
  } else if (order.method == "alternating") {
    o = rep(c(0:(p-1),(p-1):0),(max.iter+1)/2)
  } else if (order.method == "given") {
  }

  # run algorithm
  out            = caisa_order(data$X, w, sa2, pi, data$beta, r, sigma2, o,
                               max.iter, min.iter, tol$convtol, tol$epstol,
                               stepsize, update.sigma, verbose)
  
  out$intercept  = c(data$ZtZiZy - data$ZtZiZX %*% out$beta)
  out$data       = data
  
  return (out)
}