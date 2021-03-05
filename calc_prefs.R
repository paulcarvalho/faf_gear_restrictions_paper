### Author: Paul Carvalho
### Code adopted from Hall et al. (2006)
###
### Description: Calculate the lognormal probability functions for prey preferences,
###              based on prey/predator size (weight) ratio. Returns 4-dimensional
###              matrix with preferences for predator i in size class j and prey 
###              species l in size class k

calc_prefs <- function(L.mid, nsc, nspecies, mu, sigma, weight, sc_Linf){
  
  bs.ratio <- zeros(nsc, nspecies, nsc, nspecies)
  
  for(sp1 in 1:nspecies){   # predator species
    for(sp2 in 1:nspecies){ # prey species
      for(sc1 in 1:(sc_Linf[sp1])){ # only as far as the predators max size class
        for(sc2 in 1:(min(sc_Linf[sp1],sc_Linf[sp2]))){
          bs.ratio[sc1, sp1, sc2, sp2] <- (weight[sc2, sp2]) / (weight[sc1, sp1])
          bs.ratio[sc1, sp1, sc2, sp2] <- dlnorm(bs.ratio[sc1, sp1, sc2, sp2], mu, sigma)
        }
      }
    } 
  }

  return(bs.ratio)
}