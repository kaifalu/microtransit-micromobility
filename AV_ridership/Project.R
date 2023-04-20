## ----libraries, message=FALSE, warning=FALSE-----------------------------
library(reshape2)
library(plyr)
library(ggplot2)
library(xtable)
library(MASS)
library(mvtnorm)
library(BNPdensity)
library(mcmc)
library(StanHeaders)
library(rstan)
library(GGally)
library(dplyr)
library(invgamma)

## ------ Calculate log-posterior ----------------------------------
lupost_factory_weekday <- function(x, mu, sigma) {
  logd <- dlnorm(x, mu, abs(sigma), log = TRUE)
  return(sum(logd) - log(sigma^2)) # log-posterior
}

## ------ Calculate likelihood ----------------------------------
likelihood_factory_weekday <- function(x, mu, sigma) {
  logd <- dlnorm(x, mu, abs(sigma), log = TRUE)
  return(sum(logd))
}

## ----random walk metropolis hastings algorithm-------------------------
metropolis_hastings_lognorm = function(x, mu, sigma, n_points) {
  #initialize
  current = c(mu,sigma)
  mu_sigma = matrix(current, nrow=n_points, ncol=2, byrow=TRUE)
  lupost = matrix(lupost_factory_weekday(x,mu,sigma), nrow=n_points, ncol=1, byrow=TRUE)
  likelihood = matrix(likelihood_factory_weekday(x,mu,sigma), nrow=n_points, ncol=1, byrow=TRUE)
  nn = 0
  S = diag(2)
  for (i in 2:n_points) {
    proposed = mvrnorm(1,current,2.4^2*S/2)
    logr = lupost_factory_weekday(x,proposed[1],proposed[2])-lupost_factory_weekday(x,current[1],current[2])
    if (log(runif(1)) < logr) {current = proposed; nn = nn + 1}
    mu_sigma[i,] = current
    lupost[i] = lupost_factory_weekday(x,current[1],current[2])
    likelihood[i] = likelihood_factory_weekday(x,current[1],current[2])
    if (i%%50 == 0) S = var(mu_sigma[1:i,])
  }
  sampler = list(mu=mu_sigma[,1], sigma = mu_sigma[,2], log_posterior=lupost, acceptance=nn/n_points,
                 arrival_rate = exp(mu_sigma[,1]+mu_sigma[,2]^2/2),
                 ar_mean = mean(exp(mu_sigma[,1]+mu_sigma[,2]^2/2)),
                 ar_CI = c(sort(exp(mu_sigma[,1]+mu_sigma[,2]^2/2))[round(n_points*0.025)],
                           sort(exp(mu_sigma[,1]+mu_sigma[,2]^2/2))[round(n_points*0.975)]))
  return(sampler)
}

## ---------- Example of full model (intercept, X1-X4 included) ----------------

mu_MLE = 2.3629
sigma_squared_MLE = 0.7904
sigma_MLE = sqrt(sigma_squared_MLE)

## Monday
x = c(2,7,23,7,29,14,4,5,21,20,5,7,16,16,2,14,22,18,10,13,6,14,12,2,5,7,4,9,1,8,4,6,6,10,7,1)

lupost <- lupost_factory_weekday(x, mu_MLE, sigma_MLE)
print(lupost)

# sampling
model_sampler = metropolis_hastings_lognorm(x, mu_MLE, sigma_MLE, n<-1000)

print(model_sampler$acceptance)
print(model_sampler$ar_mean)
print(model_sampler$ar_CI)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2),nrow=1,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$mu, xlab="Iter", ylab=expression(alpha), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$sigma, xlab="Iter", ylab=expression(beta), main="Random Walk Metropolis-Hastings (RWMH)")

acf(model_sampler$mu,100, plot=TRUE, main = expression(mu))
acf(model_sampler$sigma,100, plot=TRUE, main = expression(sigma^2))



## Tuesday
x = c(6,4,5,19,10,21,17,6,10,3,15,9,13,4,9,8,21,4,7,15,7,4,1,9,6,9,2,47,8,32,37,6,7,3,6,13,17)

lupost <- lupost_factory_weekday(x, mu_MLE, sigma_MLE)
print(lupost)

# sampling
model_sampler = metropolis_hastings_lognorm(x, mu_MLE, sigma_MLE, n<-1000)

print(model_sampler$acceptance)
print(model_sampler$ar_mean)
print(model_sampler$ar_CI)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2),nrow=1,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$mu, xlab="Iter", ylab=expression(alpha), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$sigma, xlab="Iter", ylab=expression(beta), main="Random Walk Metropolis-Hastings (RWMH)")

acf(model_sampler$mu,100, plot=TRUE, main = expression(mu))
acf(model_sampler$sigma,100, plot=TRUE, main = expression(sigma^2))



## Wednesday
x = c(18,6,22,6,10,7,12,12,14,13,14,13,4,2,1,7,3,9,10,28,5,11,14,5,5,2,11,13,23,11,9,18,21,9,11,7,11,9,8,6)

lupost <- lupost_factory_weekday(x, mu_MLE, sigma_MLE)
print(lupost)

# sampling
model_sampler = metropolis_hastings_lognorm(x, mu_MLE, sigma_MLE, n<-1000)

print(model_sampler$acceptance)
print(model_sampler$ar_mean)
print(model_sampler$ar_CI)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2),nrow=1,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$mu, xlab="Iter", ylab=expression(alpha), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$sigma, xlab="Iter", ylab=expression(beta), main="Random Walk Metropolis-Hastings (RWMH)")

acf(model_sampler$mu,100, plot=TRUE, main = expression(mu))
acf(model_sampler$sigma,100, plot=TRUE, main = expression(sigma^2))



## Thursday
x = c(15,30,13,32,34,30,26,16,53,34,7,53,30,52,2,22,18,4,25,3,23,13,7,32,25,31,52,15,32,8,16,44,47,37,86,21,35,13,31,14,29,13,43)

lupost <- lupost_factory_weekday(x, mu_MLE, sigma_MLE)
print(lupost)

# sampling
model_sampler = metropolis_hastings_lognorm(x, mu_MLE, sigma_MLE, n<-1000)

print(model_sampler$acceptance)
print(model_sampler$ar_mean)
print(model_sampler$ar_CI)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2),nrow=1,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$mu, xlab="Iter", ylab=expression(alpha), main="Random Walk Metropolis-Hastings (RWMH)")
plot(1:n, model_sampler$sigma, xlab="Iter", ylab=expression(beta), main="Random Walk Metropolis-Hastings (RWMH)")

acf(model_sampler$mu,100, plot=TRUE, main = expression(mu))
acf(model_sampler$sigma,100, plot=TRUE, main = expression(sigma^2))
