### calc_ssb
###
### Author: Paul Carvalho
###
### Description: Calculate ssb for each function group

calc_ssb <- function(nspecies, nsc, L.lower, L.upper, N.ijt, y, Lmat){
     # calc spawning stock biomass
     ssb <- vector("numeric", length = nspecies)
     for(x in 1:nspecies){
          stock <- NULL
          for(z in 1:nsc){
               tmp_stock <- data.frame(size = seq(L.lower[z]+1,L.upper[z],by=1),
                                       stock = N.ijt[z,x,y-1]/5)
               stock <- rbind(stock, tmp_stock)
          }
          ssb[x] <- sum(subset(stock, size >= Lmat[x])$stock)
     }   
     return(ssb)
}
