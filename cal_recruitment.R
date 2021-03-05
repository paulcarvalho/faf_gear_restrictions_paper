### cal_recruitment
###
### Author: Paul Carvalho
###
### Description: Calibrate recruitment such that biomass of each functional group reaches
###              biomass/ha in remote areas of Indonesia (Campbell et al. 2020). The main
###              script calls this function with optim to find recruitment that minimizes
###              sum of squared residuals.

cal_recruitment <- function(par.in, N0, t, nsc, nspecies, M1, phi, L.lower, L.upper, W.a, W.b, Bi.F0, suit, ration, weight, sc_Linf, phi.min){
     source("calc_M2.R")
     source("calc_growth.R")
     source("calc_bio.R")
        
     r.i <- par.in[1:9]
     other <- par.in[10]
        
     N.ijt <- array(NA, dim = c(nsc, nspecies, t)) # empty shell to save data
     N.ijt[,,1] <- N0 # set initial abundance
     
     # run model
     for(i in 2:t){
          # recruitment
          N.tmp <- N.ijt[,,i-1]
          N.tmp[1,] <- N.tmp[1,] + r.i
          # calculate predation mortality
          M2 <- calc_M2(N=N.tmp, suit, ration, nspecies, nsc, other, weight, sc_Linf, phi.min)
          # natural mortality
          N.tmp <- N.tmp * exp(-(M1+M2))
          # grow
          N.tmp <- calc_growth(N.tmp, phi, nsc, nspecies)
          # NO fishing mortality for model calibration
          N.ijt[,,i] <- N.tmp
     }
     
     # calculate biomass at equilibrium
     B.ijt <- calc_bio(N.ijt, t, L.lower, L.upper, W.a, W.b, nsc, nspecies)
     B.end <- colSums(B.ijt[,,t])
     
     # calculate residual difference between model and expected
     residual <- log(B.end)-log(Bi.F0)
     
     # sum squared residuals
     res.out <- sum(residual^2)
     return(res.out)
}