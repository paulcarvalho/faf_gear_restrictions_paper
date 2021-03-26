# Functions
# Author: Paul Carvalho


calc_phi <- function(L.lower, L.upper, Linf, k, nsc){
  #' calc_phi
  #' 
  #' @description Calculates the proportion of fish from functional group i in size class j
  #' that grow to j+1 and re-scales times steps such that all fish in the fatest growing 
  #' functional group-size class comination reach the next size class. Code adapted from 
  #' Hall et al. (2006). 
  #'
  #' @param L.lower Lower limit of size class (included)
  #' @param L.upper Upper limit of size class (excluded)
  #' @param Linf Asymptotic length (cm)
  #' @param k Von Bertalanffy growth parameter
  #' @param nsc Number of size classes
  #'
  #' @return A list with (1) proportions of functional group i that grow from size class j to j+1 and (2) minimum time step
  
  # Upper and lower limits of size classes
  L.low <- repmat(L.lower, 1, length(Linf))
  L.up  <- repmat(L.upper, 1, length(Linf))
  # create matrix for VBGF parameters
  k.mat    <- repmat(t(as.matrix(k)), nsc, 1)
  Linf.mat <- repmat(t(as.matrix(Linf)), nsc, 1)
  # calculate time for fish to grow from lower to upper limit of size class
  tij           <- 1/k.mat * log((Linf.mat - L.low) / (Linf.mat - L.up)) # Hilborn and Walters 1992, p 428
  rownames(tij) <- as.character(L.lower)
  tij[which(is.infinite(tij))] <- 0 # replace Inf values with 0
  tij[which(is.na(tij))]       <- 0 # replace NaN with 0
  tij[which(tij < 0)]          <- 0 # replace negative values with 0
  
  # scale to the fastest growing species-size class combination
  phi_min <- min(tij[tij>0]) 
  for(i in 1:nspecies){
    for(j in 1:nsc){
      if(tij[j,i] > 0){
        tij[j,i] <- phi_min/tij[j,i]
      }
    }
  }
  return(list(tij, phi_min))
}

calc_growth <- function(N, phi, nsc, nspecies){
  #' calc_growth
  #'
  #' @description Calculate the number of fish that grow to next size class. Adopted from Hall et al. (2006).
  #'
  #' @param N Matrix with number of fish in size class j (row) and functional group i (col)
  #' @param phi Proportion of fish in size class j and functional group i that grow to the next size class 
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species, included in the model
  #'
  #' @return New abundance matrix (size class = row, functional group = col) after fish growth
  
  stay  <- N * (1-phi) 
  leave <- N * phi
  N.out <- stay + rbind(pracma::zeros(1, nspecies),leave[1:nsc-1,])
  return(N.out)
}

calc_ration <- function(k, Linf, nsc, nspecies, L.lower, L.upper, L.mid, W.a, W.b, phi.min, scale.Ge){
  #' calc_ration
  #' 
  #' @description Calculate the ration (ingestion rate) that must be consumed by species i in size class j 
  #' to account for modeled growth in a timestep (phi.min). Code adapted from Hall et al. (2006).
  #'
  #' @param k Von Bertalanffy growth parameter
  #' @param Linf Asymptotic length
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #' @param L.lower Lower limit of each size class (included)
  #' @param L.upper Upper limit of each size class (excluded)
  #' @param L.mid Midpoint of each size class
  #' @param W.a Length-weight conversion parameter
  #' @param W.b Length-weight conversion parameter
  #' @param phi.min Time step
  #' @param scale.Ge Scalar parameter for growth efficiency
  #'
  #' @return List with (1) ration that mught be consumed by each functional group-size class combination, (2) weight at midpoint of size class, 
  #' and (3) size class at asymptotic length (Linf)

  # Find the size class index at which each species/functional group reaches Linf
  sc_Linf <- rep(NA, nspecies)
  for(i in 1:nspecies){
    for(j in 1:nsc){
      if(Linf[i] > L.lower[j] & Linf[i] <= L.upper[j]){ sc_Linf[i] <- j }
        
    }
  }
  ration <- zeros(nrow = nsc, ncol = nspecies) # ration that must be consumed by each species/functional group-size class combination to account for growth.
  weight <- zeros(nrow = nsc, ncol = nspecies)
  g.eff  <- zeros(nrow = nsc, ncol = nspecies)  # Growth efficiency - proportion of food consumed that is converted to body mass
  
  for(i in 1:nspecies){
    for(j in 1:(sc_Linf[i]-1)){
      # Calculate the length after 1 timestep when the initial length is the mid point of a size class
      L2 <- L.mid[j] + (Linf[i] - L.mid[j]) * (1 - exp(-k[i] * phi.min))
      # Calculate growth increment from the weight at the mid-point
      weight[j,i] <- (W.a[i] * L.mid[j] ^ W.b[i])
      growth_inc  <- (W.a[i] * L2 ^ W.b[i]) - weight[j,i]
      g.eff[j,i]  <- 1 - (weight[j,i] / (W.a[i] * Linf[i] ^ W.b[i])) ^ 0.11 # equation 7 in Hall et al. (2006)
      g.eff[j,i]  <- g.eff[j,i] * 0.5 # 0.5 gives a reasonable scaling (Hall et al. 2006)
      g.eff[j,i]  <- g.eff[j,i] + (g.eff[j,i]*scale.Ge)
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
    growth_inc            <- (W.a[i] * L2 ^ W.b[i]) - weight[sc_Linf[i], i]
    g.eff[sc_Linf[i], i]  <- 1 - (weight[sc_Linf[i],i] / (W.a[i] * Linf[i] ^ W.b[i])) ^ 0.11
    g.eff[sc_Linf[i], i]  <- g.eff[sc_Linf[i], i] * 0.5
    g.eff[sc_Linf[i], i]  <- g.eff[sc_Linf[i], i] + (g.eff[sc_Linf[i], i]*scale.Ge)
    ration[sc_Linf[i], i] <- growth_inc * (1 / g.eff[sc_Linf[i], i])
  }
  output <- list(ration, weight, sc_Linf)
  return(output) # Units:grams
}

calc_prefs <- function(L.mid, nsc, nspecies, mu, sigma, weight, sc_Linf){
  #' calc_prefs
  #' 
  #' @description Calculate the lognormal probability function for prey preferecnes, based on prey/predator
  #' size (i.e., weight) ratio. Code adapted from Hall et al. (2006).
  #'
  #' @param L.mid Midpoint of each size class (length in cm)
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #' @param mu Log(mean) for prey size size preference
  #' @param sigma Log(standard deviation) for prey size preference
  #' @param weight Weight
  #' @param sc_Linf Size class at asymptotic length
  #'
  #' @return Four-dimensional array with preferences for predator i in size class j and prey species l in size class k
  bs.ratio <- zeros(nsc, nspecies, nsc, nspecies)
  
  for(sp1 in 1:nspecies){ # predator species
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

calc_suit <- function(M2_prefs, tau, nsc, nspecies, sc_Linf){
  #' calc_suit
  #'
  #' @description Calculate the suitabilities for predator of size l, functional group k 
  #' for prey of size j, functional group i. Code adapted from Hall et al. (2006).
  #'
  #' @param M2_prefs Functional group and size class specific predator preferenes for prey
  #' @param tau Foodweb matrix
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #' @param sc_Linf Size class at asymptotic length for each functional group
  #'
  #' @return Four-dimensional matrix with suitabilities for predator and prey
  
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
