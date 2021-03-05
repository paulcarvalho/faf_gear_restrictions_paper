### F_mortality
###
### Author: Paul Carvalho
###
### Description: Fishing mortality for each gear

F_mortality <- function(nsc, nspecies, fishing.effort, gear.mgmt, q){
     F.tmp <- array(NA, dim = c(nsc, nspecies, length(gear.mgmt))) # empty shell for fishing mortality of each fishing gear
     F.out <- array(0, dim = c(nsc, nspecies))
     
     # set up gear management scenario
     q.tmp <- array(NA, dim = c(nsc, nspecies, length(gear.mgmt)))
     for(i in 1:length(gear.mgmt)){
             if(gear.mgmt[i] == 1){
                     # assign catchability/selectivity and scale such that qs for all sums to one
                     q.tmp[,,i] = q[,,i]/sum(q[,,i])
             } else {
                     q.tmp[,,i] = 0
             }
     }
     
     # calculate fishing mortality
     for(i in 1:length(gear.mgmt)){
          F.tmp[,,i] <- fishing.effort[i] * q.tmp[,,i]
          F.out <- F.out + F.tmp[,,i]
     }
     
     return(F.out) # I don't think necessary to scale with time step
}


