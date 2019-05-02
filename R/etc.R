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

#phi    = posterior_quantile(fit)
#thresh = c(1e-8,0.6)
#x      = 1:length(fit$beta)
#len    = length(thresh)
#a      = posterior_quantile(fit, thresh)
#nz     = which(rowSums(a) != 0)

#df1     = data.frame(x   = as.vector(outer(thresh,rep(1,length(nz)))),
#                     y   = as.vector(t(a[nz,])),
#                     ind = as.vector(as.factor(rep(paste('x',nz, sep = ''), each = len))))
#levels(df1$ind) = paste('x',nz, sep = '')

#f1 = ggplot(df1, aes(x = x, y = y, color = ind, shape = ind)) + theme(axis.line = element_blank()) +
#  geom_line() + scale_y_continuous() +
#  geom_point() + #scale_shape_manual(values = c(1,2,3,4,5,6,8,15,16,17)) +
#  labs(x     = "posterior quantile (q)",
#       y     = "coefficients",
#       title = "The solution path of caisa")

#f1
