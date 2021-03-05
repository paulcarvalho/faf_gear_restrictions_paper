### catch
###
### Author: Paul Carvalho
###
### Description: Calculate fishing mortality

calc_catch <- function(N, F.mort){
     # calculate F mortality
     tot.catch <- N * (1 - exp(-F.mort))
     
     # Mathimatically the following won't occur, but this was added for old calculation of fishing mortliaty and kept
     if(sum(c(N - tot.catch) < 0) > 0) stop("Error: catch for one or more species-size class combinations exceeds limit")
     
     return(tot.catch)
}
     
