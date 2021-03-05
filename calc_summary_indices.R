### calc_summary_indices
###
### Author: Paul Carvalho
###
### Description: Calculate summary for model outputs

calc_summary_indices <- function(N.ijte, B.ijte, cN.ijte, cB.ijte, t, L.mid, Lmat, nspecies){
  # Calculate total biomass at equilibrium for all efforts
  B.equil <- colSums(colSums(B.ijte[,,t,]))
  
  # Calculate total catch at equilibrium for all efforts
  cB.equil <- colSums(colSums(cB.ijte[,,t,]))
  
  # Calculate biomass for each functional group at equilibrium for all efforts
  Bfg.equil <- colSums(B.ijte[,,t,])
  
  # Calculate catch for each functional group at equilibrium for all efforts
  cBfg.equil <- colSums(cB.ijte[,,t,])
  
  # Calculate mean weight. length and length at maturity
  BmeanW <- vector(mode="numeric", length=length(effort))
  BmeanL <- vector(mode="numeric", length=length(effort))
  BmeanLmat <- vector(mode="numeric", length=length(effort))
  CmeanW <- vector(mode="numeric", length=length(effort))
  CmeanL <- vector(mode="numeric", length=length(effort))
  for(e in 1:length(effort)){
    totalB <- sum(B.ijte[,,t,e])
    totalN <- sum(N.ijte[,,t,e])
    totalL <- sum(rowSums(N.ijte[,,t,e]) * L.mid)
    BmeanW[e] <- totalB/totalN
    BmeanL[e] <- totalL/totalN
    BmeanLmat[e] <- sum(colSums(N.ijte[,,t,e]) * Lmat) / totalN
    
    totalcB <- sum(cB.ijte[,,t,e])
    totalcN <- sum(cN.ijte[,,t,e])
    totalcL <- sum(rowSums(cN.ijte[,,t,e]) * L.mid)
    CmeanW[e] <- totalcB/totalcN
    CmeanL[e] <- totalcL/totalcN
  }
  
  # Calculate functional group specific mean weight, length, and length at maturity
  fg.BmeanW <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.BmeanL <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.BmeanLmat <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.CmeanW <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.CmeanL <- array(NA, dim=c(1,nspecies,length(effort)))
  for(e in 1:length(effort)){
    fg.totalB <- colSums(B.ijte[,,t,e])
    fg.totalN <- colSums(N.ijte[,,t,e])
    fg.totalL <- colSums(N.ijte[,,t,e] * repmat(as.array(L.mid),1,nspecies))
    fg.BmeanW[,,e] <- fg.totalB / fg.totalN
    fg.BmeanL[,,e] <- fg.totalL / fg.totalN
    fg.BmeanLmat[,,e] <- (colSums(N.ijte[,,t,e]) * Lmat) / fg.totalN
    
    fg.totalcB <- colSums(cB.ijte[,,t,e])
    fg.totalcN <- colSums(cN.ijte[,,t,e])
    fg.totalcL <- colSums(cN.ijte[,,t,e] * repmat(as.array(L.mid),1,nspecies))
    fg.CmeanW[,,e] <- fg.totalcB/fg.totalcN
    fg.CmeanL[,,e] <- fg.totalcL/fg.totalcN
  }
  out <- list(B.equil,cB.equil,Bfg.equil,cBfg.equil,BmeanW,BmeanL,BmeanLmat,CmeanW,CmeanL,
              fg.BmeanW,fg.BmeanL,fg.BmeanLmat,fg.CmeanW,fg.CmeanL)
}