### Author: Paul Carvalho
### Code adopted from Hall et al. (2006)
###
### Description: Calculate the predation mortality (M2) for each species in each size
###              class and return a matrix (size class x species/functional group)
###              of M2 values.

calc_M2 <- function(N, suit, ration, nspecies, nsc, other, weight, sc_Linf, phi.min){
  
  M2 <- matrix(0, nrow=nsc, ncol=nspecies)
  
  for(n in 1:nspecies){         # prey species
    for(m in 1:sc_Linf[n]){     # prey size class
      for(j in 1:nspecies){     # predator species
        for(i in 1:sc_Linf[j]){ # predator size class
          denom <- sum(suit[i,j,,] * weight * N) # denominator in equation 8 of Hall et al. (2006)
          denom <- denom + other
          if(denom > 0){
            M2[m,n] <- M2[m,n] + ration[i,j] * N[i,j] * (suit[i,j,m,n] / denom) # equation 8 of Hall et al. (2006)
          }
        }
      }
    }
  }
  
  output1 <- M2*phi.min # scale predation mortality to minimum time step (phi.min)
  
  # for(j in 1:nspecies){					# predator species
  # 	for(i in 1:sc_Linf[j]){		  # predator size class
  # 		for(n in 1:nspecies){			# prey species
  # 			for(m in 1:sc_Linf[n]){ # prey size class
  # 				denom <- sum(suit[i,j,,] * weight * N)
  # 				denom <- denom + other
  # 				if(denom > 0){
  # 					consump <- ration[i,j] * N[i,j] * (suit[i,j,m,n] * N[m,n] * weight[m,n] / denom)
  # 					M2_denom[n,j,i] <- M2_denom[n,j,i] + consump
  # 				}
  # 			}
  # 		}
  # 	}
  # }
  
  
  return(output1)
}