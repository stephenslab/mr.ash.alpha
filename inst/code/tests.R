# A short script implementing a few tests of mr.ash.alpha.
library(varbvs)

# Simulate a data set.
set.seed(1)
n          <- 400
p          <- 2000
X          <- matrix(rnorm(n*p),n,p)
beta       <- rep(0,p)
beta[1:10] <- rnorm(10)
y          <- drop(-2 + X %*% beta + rnorm(n))

# Try the different orderings.
fit <- glmnet::glmnet(X,y)
i1  <- univar.order(X,y)
i2  <- path.order(fit)
i3  <- absolute.order(coef(fit)[-1])

# Fit the mr.ash model, and compute posterior expectations of
# interest (means, variances, and posterior assignment probabilities).
fit1 <- mr.ash(X,y)
out  <- get.full.posterior(fit1)

# Compare against varbvsmix solution.
fit2 <- varbvsmix(X,NULL,y,alpha = matrix(0.5,p,20),mu = matrix(0,p,20),
                  sa = fit1$data$sa2,update.sigma = TRUE,
                  update.sa = FALSE,update.w = TRUE,verbose = FALSE)
plot(fit1$beta,rowSums(fit2$alpha * fit2$mu),pch = 20)
plot(fit1$pi + 1e-4,fit2$w + 1e-4,pch = 20,log = "xy")

# Predict the regression outcomes in the training data, and compare
# against the dround-truth values.
yest <- predict(fit1,X)
plot(y,yest,pch = 20,col = "darkblue")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
