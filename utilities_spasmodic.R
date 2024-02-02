# Functions to determine if a stock is spasmodic
# Original code by Paula Silvar Villamidou



## empirical cumulative distribution function
## could use R's ecdf function but to extract values
## requires defining and applying the function

## own ecdf function
ecdf_fn <- function(x){
  sx <- sort(x)
  n <- length(x)
  return(list(x = sx, y = (1:n)/n))
}


## simulate pointwise quantiles of lognormal cdf with given variability
get_bounds <- function(n, sd, alpha = 0.2, m){
  ##-------------------------------------
  ## simulates cdfs from a scaled lognormal distribution
  ## n is the length of the recruitment timeseries
  ## sd is the standard deviation on the log scale of a lognormal distributions
  ## alpha is the significance level of the bands
  ## m is the number of replicates
  ##-------------------------------------
  ## set up a container for the results
  res <- data.frame(x = rep(NA, m*n), y = rep(NA, m*n))
  ## simulate scaled cdfs m times and store
  for(j in 1:m){
    x <- rlnorm(n, meanlog = 5, sdlog = sd) ## mean doesn't matter here
    x <- x / max(x)
    ## ecdf
    ecdfj <- ecdf_fn(x)
    ## store the simulated ecdf
    res$y[((j-1)*n + 1):(j*n)] <- ecdfj$y
    ## store the scaled x
    ## rounding here to aggregate subsequently
    ## could change precision and increase m for smoother bounds
    res$x[((j-1)*n + 1):(j*n)] <- round(ecdfj$x, 2)
  }

  
  ## get the bounds
  ## lower
  lwr <- aggregate(y ~ x, quantile, p = alpha/2, data = res)
  names(lwr)[names(lwr) == "y"] <- "lwr"
  ## upper
  upr <- aggregate(y ~ x, quantile, p = 1 - alpha/2, data = res)
  names(upr)[names(upr) == "y"] <- "upr"
  bounds <- merge(lwr, upr)
  ##
  return(bounds)
}