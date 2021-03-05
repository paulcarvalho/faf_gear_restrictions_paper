### calc_sensitivity_indices
###
### Author: Paul Carvalho
###
### Description: Calculate indicators for sensitivity to model parameters

calc_sensitivity_indices <- function(N.ijte, B.ijte, cN.ijte, cB.ijte, L.mid, Lmat, tau, base.effort){
  # Calculate total biomass at equilibrium for all efforts
  B.equil <- colSums(colSums(B.ijte[,,t,]))
  
  # Calculate total catch at equilibrium for all efforts
  cB.equil <- colSums(colSums(cB.ijte[,,t,]))
  
  # Calculate mean length of stock and catch
  BmeanL <- vector(mode="numeric", length=length(effort))
  CmeanL <- vector(mode="numeric", length=length(effort))
  for(e in 1:length(base.effort)){
    totalN <- sum(N.ijte[,,t,e])
    totalL <- sum(rowSums(N.ijte[,,t,e]) * L.mid)
    BmeanL[e] <- totalL/totalN
    
    totalcN <- sum(cN.ijte[,,t,e])
    totalcL <- sum(rowSums(cN.ijte[,,t,e]) * L.mid)
    CmeanL[e] <- totalcL/totalcN
  }
  
  # Create dataframe
  si.df <- data.frame(effort = base.effort,
                      B = B.equil,
                      cB = cB.equil,
                      BmeanL = BmeanL,
                      CmeanL = CmeanL)
  
  # Find values at 0.5B0
  si.df <- si.df %>%
    mutate(sq.diff = (B - (max(B)/2))^2)
  out <- si.df[which(si.df$sq.diff == min(si.df$sq.diff)),]
  
  return(out)
}
