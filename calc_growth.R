### calc_growth
###
### Author: Paul Carvalho
###
### Description: Calculate the number of fish that grow to next size class.
###              Adopted from Hall et al. (2006)

calc_growth <- function(N, phi, nsc, nspecies){
     stay <- N * (1-phi) 
     leave <- N * phi
     N.out <- stay + rbind(pracma::zeros(1, nspecies),leave[1:nsc-1,])
     return(N.out)
}