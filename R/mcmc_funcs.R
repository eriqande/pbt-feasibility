

#' simulate a dirichlet R.V.
#' 
#' @param a  The parameter vector
rdirichlet <- function(a) {
  y <- rgamma(n = length(a), shape = a)
  y / sum(y)
}