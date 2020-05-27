#' @rdname mr.ash
#' 
#' @export
#' 
mr.ash.dev                  = function(X, y, Z = NULL, sa2 = NULL,
                                       method = c("caisa","accelerate","block",
                                                  "sigma","sigma_scaled",
                                                  "sigma_indep"),
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
  if (is.null(sa2)) {
    
  } else {
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
    } else if (method == "sigma_scaled") {
      out          = caisa_em2  (data$y, data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
      out$beta     = out$beta * sqrt(out$sigma2)
    } else if (method == "sigma_indep") {
      out          = caisa_em3  (data$X, w, sa2, pi, data$beta, r, sigma2,
                                 max.iter, min.iter, tol$convtol, tol$epstol,
                                 update.sigma2, verbose)
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
    warning(paste0("the estimated mixtuer proportion pi[",K,"] is non-zero -- consider increasing the range of grid variances until sa2[",K,"] is zero\n",sep = ""))
  }
  
  return(out)
}
