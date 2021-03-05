### Author: Paul Carvalho
### Code adopted from Hall et al. (2006)
###
### Description: Calculate the ration (ingestion rate) that must be consumed
###              by species i in size class j to account for modeled growth
###              in a timestep (phi.min)

calc_ration <- function(k, Linf, nsc, nspecies, L.lower, L.upper, L.mid, W.a, W.b, phi.min, scale.Ge){
  
  # Find the size class index at which each species/functional group reaches Linf
  sc_Linf <- rep(NA, nspecies)
  for(i in 1:nspecies){
    for(j in 1:nsc){
      
      if(Linf[i] > L.lower[j] & Linf[i] <= L.upper[j])
        sc_Linf[i] = j
    }
  }
  
  ration <- zeros(nrow = nsc, ncol = nspecies) # ration that must be consumed by each species/functional group-size class combination to account for growth.
  weight <- zeros(nrow = nsc, ncol = nspecies)
  g.eff <- zeros(nrow = nsc, ncol = nspecies)  # Growth efficiency - proportion of food consumed that is converted to body mass
  
  for(i in 1:nspecies){
    for(j in 1:(sc_Linf[i]-1)){
      
      # Calculate the length after 1 timestep when the initial length is the mid point of a size class
      L2 <- L.mid[j] + (Linf[i] - L.mid[j]) * (1 - exp(-k[i] * phi.min))
      
      # Calculate growth increment from the weight at the mid-point
      weight[j,i] <- (W.a[i] * L.mid[j] ^ W.b[i])
      growth_inc <- (W.a[i] * L2 ^ W.b[i]) - weight[j,i]
      g.eff[j,i] <- 1 - (weight[j,i] / (W.a[i] * Linf[i] ^ W.b[i])) ^ 0.11 # equation 7 in Hall et al. (2006)
      g.eff[j,i] <- g.eff[j,i] * 0.5 # 0.5 gives a reasonable scaling (Hall et al. 2006)
      g.eff[j,i] <- g.eff[j,i] + (g.eff[j,i]*scale.Ge)
      ration[j,i] <- growth_inc * (1 / g.eff[j,i]) # ingestion rate that must be consumed by species/functional group i in size class j to account for modeled growth in a given timestep (tonnes/timestep)
    }
  }
  
  # repeat above for largest size class
  for(i in 1:nspecies){
    # Calculate the mid-point between the lower bound of the largest size class and the Linf
    fmid <- L.lower[sc_Linf[i]] + ((Linf[i] - L.lower[sc_Linf[i]]) / 2)
    
    # Calculate the length after 1 timestep when the initial length is the mid-point of a size class
    L2 <- fmid + (Linf[i] - fmid) * (1 - exp(-k[i] * phi.min))
    
    # Calculate growth increment from the weight at the mid-point
    weight[sc_Linf[i], i] <- W.a[i] * fmid ^ W.b[i]
    growth_inc = (W.a[i] * L2 ^ W.b[i]) - weight[sc_Linf[i], i]
    g.eff[sc_Linf[i], i] <- 1 - (weight[sc_Linf[i],i] / (W.a[i] * Linf[i] ^ W.b[i])) ^ 0.11
    g.eff[sc_Linf[i], i] <- g.eff[sc_Linf[i], i] * 0.5
    g.eff[sc_Linf[i], i] <- g.eff[sc_Linf[i], i] + (g.eff[sc_Linf[i], i]*scale.Ge)
    ration[sc_Linf[i], i] <- growth_inc * (1 / g.eff[sc_Linf[i], i])
  }
  
  output <- list(ration, weight, sc_Linf)
  return(output) # Units:grams
}
