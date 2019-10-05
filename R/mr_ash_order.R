#' @title Something
#' 
#' @description
#' 
#' @details Add details here.
#'
#' @param X Describe X here. Sparse matrix? Dense?
#' 
#' @param y Describe y here.
#' 
#' @return Describe return value here.
#' 
#' @useDynLib varbvs2
#' 
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' 
#' @examples
#' 
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y = X %*% beta + rnorm(n)
#' 
#' fit.mr_ash  = mr_ash_order(X,y, order.method = "random")
#' 
#' # DO NOT RUN
#' 
#' # fit.mr_ash  = mr_ash_order(X,y, order.method = "normal")
#' 
#' # fit.mr_ash  = mr_ash_order(X,y, order.method = "given", o = rep(0:(p-1), max.iter))
#' 
#' Xnew        = matrix(rnorm(n*p),n,p)
#' ynew        = predict.mr_ash(fit.mr_ash, Xnew)
#' 
#' 
#' @export
#' 
mr_ash_order             = function(X, y, Z = NULL, sa2 = NULL,
                                    order.method = "increasing", o = NULL, K = 10,
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
  
  if (order.method == "increasing") {
    o = rep(1:p, max.iter)
  } else if (order.method == "decreasing") {
    o = rep(p:1, max.iter)
  } else if (order.method == "random") {
    o = sample(1:p, p)
    for (i in 2:max.iter) {
      o = c(o, sample(1:p, p))
    }
  } else if (order.method == "alternating") {
    o = rep(c(1:p,p:1),(max.iter+1)/2)
  } else if (order.method == "manual") {
  }
  update.order = o - 1;

  # run algorithm
  out            = caisa_order(data$X, w, sa2, pi, data$beta, r, sigma2, update.order,
                               max.iter, min.iter, tol$convtol, tol$epstol,
                               stepsize, update.sigma, verbose)
  
  out$intercept  = c(data$ZtZiZy - data$ZtZiZX %*% out$beta)
  out$data       = data
  out$order      = update.order
  
  class(out)      <- c("mr_ash","list")
  return (out)
}