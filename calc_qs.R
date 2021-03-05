### Author: Paul Carvalho
###
### Description: Calculate the product of size selectivity and
###			  functional group catchability
calc_qs <- function(x, mu, sd, q){
  qs <- (as.matrix(dlnorm(x, mu, sd))) %*% t(as.matrix(q))
  qs[qs < 0] <- 0 
  qs <- qs / sum(qs)
  return(qs)
}