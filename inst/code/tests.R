# A short script implementing a few tests of mr.ash.alpha.

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

# Predict the regression outcomes in the training data, and compare
# against the dround-truth values.
yest <- predict(fit1,X)
plot(y,yest,pch = 20,col = "darkblue")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
