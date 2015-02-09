


#' simulate a compound dirichlet multinomial
#' 
#' The a parameter recycles so it is easy to simulate one with constant a 
#' 
#' @param a The alpha parameter
#' @param K Number of categories
#' @param size sample size
#' @param n size number of replicates to simulate
simCDM <- function(n, size, a, K) {
  a <- rep(a, K)
  gams <- rgamma(K, shape = a, scale = 100)
  p = gams/sum(gams)
  
  rmultinom(n, size, p)
}


#' given a matrix of simCDM reps, return the sd and the number that are non-zero
#' 
cdm_stats <- function(M) {
  # note that we can just average over all of them
  ret <- vector()
  ret["sd"] <- sd(M[M>0])
  ret["nz"] <- sum(M>0) / ncol(M)
  ret
}


#' here is a function that returns a "score" for closeness to sd and non-zero
#' 
#' @param x is a vector. component 1 is sd and 2 is nz
#' @param p is a vector like x, but is the fixed value we are trying to get close to
dist_cdm <- function(x, p) {
  sqrt(sum(((x - p) / p) ^ 2))
}


#' and here is an objective function given a and K as elements of the first arg
obj <- function(val, true_sd, true_nz, size, reps = 10) {
  x <- cdm_stats(
    simCDM(reps, size = size, a = val[1], K = val[2])
  )
  
  dist_cdm(x, c(true_sd, true_nz))
}