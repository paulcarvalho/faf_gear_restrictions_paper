### calc_bio
###
### Author: Paul Carvalho
###
### Description: Convert length to biomass using the midpoint of each size class

calc_bio <- function(N, t, L.lower, L.upper, W.a, W.b, nsc, nspecies){
     
     # Save N as an array if it's a matrix
     if(t == 1){ 
        N <- array(N, dim=c(nsc,nspecies,t))
     }
        
     # calculate midpoint of each size class
     Lmid <- (L.lower + L.upper) / 2
     Lmid.mat <- pracma::repmat(as.matrix(Lmid), 1, nspecies)     
     Lmid.array <- array(Lmid.mat, dim=c(nsc,nspecies,t))
     
     # create matrices for params
     wa.mat <- pracma::repmat(t(as.matrix(W.a)), nsc, 1)
     wa.array <- array(wa.mat, dim=c(nsc,nspecies,t))
     wb.mat <- pracma::repmat(t(as.matrix(W.b)), nsc, 1)
     wb.array <- array(wb.mat, dim=c(nsc,nspecies,t))
     
     # calculate biomass (kg)
     N.out <- N * ((wa.array*Lmid.array^wb.array) / 1000)
     return(N.out)
}

