#' @title Mr.ASH (Multiple Regression with Adaptive Shrikage)
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
#' @useDynLib mr.ash.alpha
#' 
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' @importFrom stats var
#' 
#' @examples
#' 
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' fit.mr.ash  = mr.ash(X,y, method = "caisa")
#' # fit.mr.ash  = mr.ash(X,y, method = "update_g")
#' # fit.mr.ash  = mr.ash(X,y, method = "accelerate")
#' 
#' Xnew        = matrix(rnorm(n*p),n,p)
#' ynew        = Xnew %*% beta + rnorm(n)
#' ypred       = predict(fit.mr.ash, Xnew)

#' rmse        = norm(ynew - ypred, '2') / sqrt(n)
#' betahat     = fit.mr.ash$beta
#' 
#' @export
#' 
mr.ash                      = function(X, y, Z = NULL, sa2 = c((2^((0:19)/20) - 1)^2, NULL),
                                       K = 20, method = c("caisa","accelerate","block","sigma_scaled","sigma_indep"),
                                       max.iter = 1000, min.iter = 1,
                                       beta.init = NULL,
                                       update.pi = TRUE, pi = NULL,
                                       update.sigma = TRUE, sigma2 = NULL,
                                       update.order = NULL,
                                       standardize = FALSE, intercept = TRUE,
                                       mixsqpiter = 5, mode = FALSE,
                                       tol = set_default_tolerance(),
                                       verbose = TRUE){
  
  # get sizes
  n            = dim(X)[1];
  p            = dim(X)[2];
  
  # match method
  method      <- match.arg(method)
  
  # set default tolerances unless specified
  tol0         = set_default_tolerance()
  tol          = modifyList(tol0,tol,keep.null = TRUE)
  
  # calculate sa2 if missing
  if ( is.null(sa2) ) {
    betahat    = (t(X) %*% y) / n
    sa2        = logspace( b = 2 * ceiling(max(betahat^2)), length = K)
  }
  K            = length(sa2);
  
  # remove covariates
  data         = remove_covariate(X, y, Z, standardize, intercept);
  data$sa2     = sa2
  
  # initialize beta
  if ( is.null(beta.init) ){
    data$beta    = as.vector(double(p));
    # initialize r
    r            = data$y
  } else {
    # scale if needed
    if (standardize) {
      data$beta  = as.vector(beta.init) * attr(data$X,"scaled:scale");
    } else {
      data$beta  = as.vector(beta.init);
    }
    # initialize r
    r            = data$y - data$X %*% data$beta;
  }
  data$beta[1] = data$beta[1] + 0;      # to make sure beta.init is not modified
  
  # sigma2
  if ( is.null(sigma2) ) {
    sigma2 = c(var(r));
  }
  
  # precalculate
  w            = colSums(data$X^2);
  data$w       = w;
  
  # initialize other parameters
  if ( is.null(pi) ) {
    if ( is.null(beta.init) ){
      
      Phi          = matrix(1,p,K)/K;
      pi           = rep(1,K)/K;
      
    } else {
      
      S            = outer(1/w, sa2, '+') * sigma2;
      Phi          = -data$beta^2/S/2 - log(S)/2;
      Phi          = exp(Phi - apply(Phi,1,max));
      Phi          = Phi / rowSums(Phi);
      pi           = colMeans(Phi);
      
    }
  } else {
    Phi          = matrix(rep(pi, each = p), nrow = p)
  }
  
  # run algorithm
  
  if ( is.null(update.order) ) {
    
    update.order   = 1:p
    if (method == "caisa") {
      if (update.pi) {
        out          = caisa_em   (data$X, w, sa2, pi, data$beta, r, sigma2,
                                   max.iter, min.iter, tol$convtol, tol$epstol,
                                   update.sigma, verbose)
      } else {
        out          = caisa_fix_pi(data$X, w, sa2, pi, data$beta, r, sigma2,
                                    max.iter, min.iter, tol$convtol, tol$epstol,
                                    update.sigma, verbose)
      }
    } else if (method == "accelerate") {
      out           = caisa_acc (data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, mixsqpiter, tol$convtol, tol$epstol,
                                 update.sigma, verbose)
    } else if (method == "block") {
      stepsize      = 1
      out           = caisa_g  (data$X, w, sa2, Phi, pi, data$beta, r, sigma2,
                                max.iter, min.iter, tol$convtol, tol$epstol,
                                stepsize, update.sigma, mode, verbose)
    } else if (method == "sigma_scaled") {
      out          = caisa_em2  (data$y, data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma, verbose)
      out$beta     = out$beta * sqrt(out$sigma2)
    } else if (method == "sigma_indep") {
      out          = caisa_em3  (data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma, verbose)
    }
<<<<<<< HEAD
  } else if (update.order == "random") { # order for random permutation
    update.order    = sample(0:(p-1))
    for (i in 2:max.iter) {
      update.order  = c(update.order, sample(0:(p-1)))
    }
=======
  } else if (is.numeric(update.order)) {
    update.order    = rep(update.order - 1, max.iter)
>>>>>>> 8776f73558cd8fc94653965259fdef2398a579a9
    out             = caisa_order(data$X, w, sa2, pi, data$beta, r, sigma2, update.order,
                                  max.iter, min.iter, tol$convtol, tol$epstol,
                                  update.sigma, verbose)
  } else if (update.order == "random") { # order for random permutation
<<<<<<< HEAD
    update.order= sample(0:(p-1))
    for (i in 2:max.iter) {
      update.order  = c(update.order, sample(0:(p-1)))
=======
    update.order    = sample(0:(p-1))
    for (i in 2:max.iter) {
      update.order = c(update.order, sample(0:(p-1)))
>>>>>>> 8776f73558cd8fc94653965259fdef2398a579a9
    }
    out             = caisa_order(data$X, w, sa2, pi, data$beta, r, sigma2, update.order,
                                  max.iter, min.iter, tol$convtol, tol$epstol,
                                  update.sigma, verbose)
  }
  
  # return intercept, processed data and the update order
  out$intercept    = c(data$ZtZiZy - data$ZtZiZX %*% out$beta)
  out$data         = data
  out$update.order = update.order
  
  # if standardize, then rescale beta
  if (standardize) {
    out$beta       = out$beta / attr(data$X,"scaled:scale")
    
  }
  out$pi           = drop(out$pi)
  out$beta         = drop(out$beta)
  
  class(out)      <- c("mr.ash","list")
  return (out)
}
