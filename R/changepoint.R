# this code is heavily imported from Matthew Stephens susier investigation for changepoint problem
# https://stephens999.github.io/susier-investigate/changepoint.html

# wrapper for caisa
caisa_cp = function(y, X, sa2 = NULL, init = NULL, ...){

  out = varmixopt(X = X, Z = NULL, y = y, stepsize = 1, max.iter = 2000, sa2 = sa2, method = "update_g",
                  update.sigma = TRUE, beta.init = init, tol = list(epstol = 1e-14, convtol = 1e-8, ...))
  
  yhat = predict.caisa(out, X)
  
  return(list(fit = out, yhat = yhat))
}

# wrapper for bcp
bcp_cp = function(y) {
  fit.bcp = bcp(y)
  yhat.bcp = fit.bcp$posterior.mean
  
  return(list(fit = fit.bcp, yhat = yhat.bcp))
}

# wrapper for bayesb
bayesb_cp = function(y, X) {
  fit.bayesb = BGLR(y, ETA = list(list(X = X, model = "BayesB")), verbose = FALSE)
  yhat.bayesb = fit.bayesb$yHat
  
  return(list(fit = fit.bayesb, yhat = yhat.bayesb))
}

# wrapper for susier
susie_cp = function(y,X,auto=FALSE,standardize=FALSE,...){

  if(auto){
    s = susie_auto(X,y,standardize=standardize,...)
  } else {
    s = susie(X,y,standardize=standardize,...)
  }
  yhat.s = predict(s)
  return(list(fit = s, yhat = yhat.s))
}

# wrapper for L0Learn
l0_cp = function(y,algorithm="CDPSI",maxSuppSize=20,...){ 
  n=length(y)
  X = matrix(0,nrow=n,ncol=n-1)
  for(j in 1:(n-1)){
    for(i in (j+1):n){
      X[i,j] = 1
    }
  }
  y.l0.cv = L0Learn.cvfit(X,y,nFolds=5,seed=1,penalty="L0",maxSuppSize = maxSuppSize,algorithm=algorithm,...) 
  opt = which.min(y.l0.cv$cvMeans[[1]])
  yhat = predict(y.l0.cv, newx = X,lambda=y.l0.cv$fit$lambda[[1]][opt])@x
  return(list(fit = y.l0.cv$fit,yhat=yhat))
}

# wrapper for trendfiltering
tf_cp = function(x){
  x.tf = trendfilter(x,ord=0)
  x.tf.cv = cv.trendfilter(x.tf)
  opt = which(x.tf$lambda==x.tf.cv$lambda.min) #optimal value of lambda
  yhat= x.tf$fit[,opt]
  return(list(fit=x.tf,yhat =yhat))
}

# wrapper for 
segment_cp = function(x){
  res = segment(CNA(x,rep(1,length(x)),1:length(x)))
  yhat = rep(res$output$seg.mean,diff(c(0,res$output$loc.end)))
  return(list(fit=res,yhat=yhat))
}