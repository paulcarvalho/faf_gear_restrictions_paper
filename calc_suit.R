### Author: Paul Carvalho
### Code adopted from Hall et al. (2006)
###
### Description: Calculate the suitabilities for predator of size l,
###              species k of prey of size j, species i.

calc_suit <- function(M2_prefs, tau, nsc, nspecies, sc_Linf){
  
  suit <- zeros(nsc, nspecies, nsc, nspecies)
  
  for(i in 1:nspecies){ # prey species
    for(j in 1:nsc){    # prey size class
      for(l in 1:nsc){  # predator size class
        suit[l,,j,i] <- M2_prefs[l,,j,i] * tau[,i] # equation 8.1 in Hall et al. (2006)
      }
    }
  }
  
  # standardize such that suitabilities sum to one
  # if the species is not a predator of anything, denom is zero
  denom <- matrix(NA, nsc, nspecies)
  for(sp1 in 1:nspecies){
    for(sc1 in 1:nsc){
      denom[sc1, sp1] <- sum(sum(suit[sc1,sp1,,]))
      if(denom[sc1, sp1] > 0){
        suit[sc1,sp1,,] <- suit[sc1,sp1,,] / denom[sc1, sp1]
      }
    }
  }
  
  return(suit)
}