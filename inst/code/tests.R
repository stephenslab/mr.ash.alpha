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
out <- univar.order(X,y)

# Fit the mr.ash model.
fit1 <- mr.ash(X,y)

# Predict the regression outcomes in the training data, and compare
# against the dround-truth values.
yest <- predict(fit1,X)
plot(y,yest,pch = 20,col = "darkblue")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
