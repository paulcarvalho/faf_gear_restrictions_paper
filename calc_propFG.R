### calc_propFG
###
### Author: Paul Carvalho
###
### Description: Calculate the proportion of biomass in each functional group

calc_propFG <- function(N_bio){
     sum_bio <- colSums(N_bio)
     prop_bio <- sum_bio / sum(sum_bio)
     return(prop_bio)
}