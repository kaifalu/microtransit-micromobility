## ----libraries, message=FALSE, warning=FALSE-----------------------------
library(plyr)
library(ggplot2)
library(xtable)
library(MASS)
library(mvtnorm)
library(mcmc)
library(GGally)
library(dplyr)

## ------ Calculate log-posterior ----------------------------------
lupost_factory_weekend <- function(x, df, scale) {
  df=abs(df)
  scale = abs(scale)
  logd <- -x^2/(2*scale^2)+(df-1)*log(x)-(df/2-1)*log(2)-(df-1)*log(scale)-log(gamma(df/2))
  return(sum(logd) - log(1+scale^2)) # log-posterior
}

## ------ Calculate likelihood ----------------------------------
likelihood_factory_weekend <- function(x, df, scale) {
  df=abs(df)
  scale = abs(scale)
  logd <- -x^2/(2*scale^2)+(df-1)*log(x)-(df/2-1)*log(2)-(df-1)*log(scale)-log(gamma(df/2))
  return(sum(logd))
}

## ----random walk metropolis hastings algorithm-------------------------
metropolis_hastings_chi = function(x, df, scale, n_points) {
  #initialize
  current = df
  df_scale = matrix(current, nrow=n_points, ncol=1, byrow=TRUE)
  lupost = matrix(lupost_factory_weekend(x,df,scale), nrow=n_points, ncol=1, byrow=TRUE)
  likelihood = matrix(likelihood_factory_weekend(x,df,scale), nrow=n_points, ncol=1, byrow=TRUE)
  nn = 0
  S = 1
  for (i in 2:n_points) {
    proposed = rnorm(1,current,S)
    proposed = abs(proposed)
    logr = lupost_factory_weekend(x,proposed[1],scale)-lupost_factory_weekend(x,current[1],scale)
    if (log(runif(1)) < logr) {current = proposed; nn = nn + 1}
    df_scale[i] = current
    lupost[i] = lupost_factory_weekend(x,current[1],scale)
    likelihood[i] = likelihood_factory_weekend(x,current[1],scale)
    if (i%%50 == 0) S = var(df_scale[1:i])
  }
  sampler = list(df=df_scale, log_posterior=lupost, acceptance=nn/n_points,
                 arrival_rate = sqrt(2)*scale*gamma((df_scale+1)/2)/gamma(df_scale/2),
                 ar_mean = mean(sqrt(2)*scale*gamma((df_scale+1)/2)/gamma(df_scale/2)),
                 ar_CI = c(sort(sqrt(2)*scale*gamma((df_scale+1)/2)/gamma(df_scale/2))[round(n_points*0.025)],
                           sort(sqrt(2)*scale*gamma((df_scale+1)/2)/gamma(df_scale/2))[round(n_points*0.975)]))
  return(sampler)
}

## ---------- Example of full model (intercept, X1-X4 included) ----------------

df_MLE = 1.9907
scale_MLE = 55.3774

## Friday
x = c(68,13,64,62,56,51,17,60,133,43,208,20,44,38,47,4,60,41,41,21,25,20,35,28,11,57,57,55,7,37,34,55,61,35,56,95,69,42,44,25,79,46,61,47,72)

lupost <- lupost_factory_weekend(x, df_MLE, scale_MLE)
print(lupost)

# sampling
model_sampler = metropolis_hastings_chi(x, df_MLE, scale_MLE, n<-1000)

print(model_sampler$acceptance)
print(model_sampler$ar_mean)
print(model_sampler$ar_CI)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2),nrow=1,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$df, xlab="Iter", ylab=expression(k), main="Friday: RWMH")
#plot(1:n, model_sampler$scale, xlab="Iter", ylab=expression(s), main="Random Walk Metropolis-Hastings (RWMH)")

#acf(model_sampler$df,100, plot=TRUE, main = expression(k))
#acf(model_sampler$scale,100, plot=TRUE, main = expression(s))



## Satursday
x = c(135,108,86,127,153,151,99,99,105,117,44,73,66,109,100,87,90,84,67,76,43,71,149,96,63,89,58,90,48,96,57,107,123,129,114,114,61,37,37,92,134,85,60)

lupost <- lupost_factory_weekend(x, df_MLE, scale_MLE)
print(lupost)

# sampling
model_sampler = metropolis_hastings_chi(x, df_MLE, scale_MLE, n<-1000)

print(model_sampler$acceptance)
print(model_sampler$ar_mean)
print(model_sampler$ar_CI)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2),nrow=1,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$df, xlab="Iter", ylab=expression(k), main="Saturday: RWMH")
#plot(1:n, model_sampler$scale, xlab="Iter", ylab=expression(s), main="Random Walk Metropolis-Hastings (RWMH)")

#acf(model_sampler$df,100, plot=TRUE, main = expression(k))
#acf(model_sampler$scale,100, plot=TRUE, main = expression(s))



## Sunday
x = c(87,87,90,120,83,97,90,23,77,71,46,69,21,46,33,53,164,16,77,64,27,25,69,53,52,102,29,97,87,97,96,45,87,77,11,135,67,76,39,55,63,70,20,91,54,64)

lupost <- lupost_factory_weekend(x, df_MLE, scale_MLE)
print(lupost)

# sampling
model_sampler = metropolis_hastings_chi(x, df_MLE, scale_MLE, n<-1000)

print(model_sampler$acceptance)
print(model_sampler$ar_mean)
print(model_sampler$ar_CI)

# plot the coefficients by iteration
F_m1 <- matrix(c(1,2),nrow=1,ncol=2,byrow=TRUE)
layout(F_m1)
plot(1:n, model_sampler$df, xlab="Iter", ylab=expression(k), main="Sunday: RWMH")
#plot(1:n, model_sampler$scale, xlab="Iter", ylab=expression(s), main="Random Walk Metropolis-Hastings (RWMH)")

#acf(model_sampler$df,100, plot=TRUE, main = expression(k))
#acf(model_sampler$scale,100, plot=TRUE, main = expression(s))
