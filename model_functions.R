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
  L.up <- repmat(L.upper, 1, length(Linf))
  
  # create matrix for VBGF parameters
  k.mat <- repmat(t(as.matrix(k)), nsc, 1)
  Linf.mat <- repmat(t(as.matrix(Linf)), nsc, 1)
  
  # calculate time for fish to grow from lower to upper limit of size class
  tij <- 1/k.mat * log((Linf.mat - L.low) / (Linf.mat - L.up)) # Hilborn and Walters 1992, p 428
  rownames(tij) <- as.character(L.lower)
  tij[which(is.infinite(tij))] <- 0 # replace Inf values with 0
  tij[which(is.na(tij))] <- 0       # replace NaN with 0
  tij[which(tij < 0)] <- 0          # replace negative values with 0
  
  # scale to the fastest growing species-size class combination
  phi_min <- min(tij[tij>0]) 
  for(i in 1:nspecies){
    for(j in 1:nsc){
      if(tij[j,i] > 0){
        tij[j,i] = phi_min/tij[j,i]
      }
    }
  }
  return(list(tij, phi_min))
}
