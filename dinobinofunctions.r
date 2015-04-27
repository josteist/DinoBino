## Functions:
# To find the mode of a list of numbers
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Drawing a number/index from an empirical probability distribution.
Emprand <- function(x) {
  j <- runif(1) ## drawing random uniform number
  which(cumsum(x)/sum(x)>j)[1] #returning first entry in the cumsum/sum [empirical density function]
}

## Likelihood function
likepoissint <- function(lambda,dt,nobs) {
  (((lambda*dt)^(nobs))/(factorial(nobs))*(exp(-lambda*dt)))/(1- exp(-lambda*dt))
}

## Negative Log likelihood for many observations
loglikepoissint <- function(lambda,dt,nobs) {
  ll <- 0
  for (ii in (1:length(dt))){
    ll <- ll - log(likepoissint(lambda,dt[ii],nobs[ii]))
  }
  return(ll)
}
