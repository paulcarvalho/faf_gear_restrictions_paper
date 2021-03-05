### catchability
###
### Author: Paul Carvalho
###
### Description: Calculates the catchability of functional group-size class
###              combinations based on Manly-Chesson selectivity indices or
###              probabilities from logistic regression model
###              (Carvalho et al. in prep).

catchability <- function(nsc, nspecies, q.option = 1, ci = 1){
     if(q.option == 1){
          # load data
          sel_fg <- read_xlsx("model_params.xlsx", sheet = 2)
          sel_size <- read_xlsx("model_params.xlsx", sheet = 3)
          
          # convert to matrix form
          fg <- as.matrix(sel_fg[,2:4])
          size <- as.matrix(sel_size[,2:4])
          q <- array(NA, dim = c(nsc, nspecies, ncol(fg))) # D1(rows) = size classes, D2(cols) = func groups, D3 = fishing gears
                                                           # D3: 1-line, 2-net, 3-spear
          
          # calculate catchability
          for(i in 1:ncol(fg)){
               gamma_fg <- t(as.matrix(fg[,i]))
               gamma_size <- as.matrix(size[,i])
               q[,,i] <- gamma_size %*% gamma_fg
               q[,,i] <- q[,,i] / sum(q[,,i])
          }
          
     } else if (q.option == 2) {
          # load data
          sel_fg <- read_xlsx("model_params.xlsx", sheet = 4)
          sel_size <- read_xlsx("model_params.xlsx", sheet = 5)
          
          if(ci == 1){
               fg <- as.matrix(round(sel_fg[,c(2,5,8)],3))
               size <- as.matrix(round(sel_size[,c(2,5,8)],3))
          } else if(ci == 2){
               fg <- as.matrix(round(sel_fg[,c(3,6,9)],3))
               size <- as.matrix(round(sel_size[,c(3,6,9)],3))
          } else if(ci == 3){
               fg <- as.matrix(round(sel_fg[,c(4,7,10)],3))
               size <- as.matrix(round(sel_size[,c(4,7,10)],3))
          }
          q <- array(NA, dim = c(nsc, nspecies, ncol(fg))) # D1(rows) = size classes, D2(cols) = func groups, D3 = fishing gears
                                                           # D3: 1-line, 2-net, 3-spear
          
          # calculate catchability
          for(i in 1:ncol(fg)){
               gamma_fg <- t(as.matrix(fg[,i]))
               gamma_size <- as.matrix(size[,i])
               q[,,i] <- gamma_size %*% gamma_fg
               q[,,i] <- q[,,i] / sum(q[,,i])
          }
     }
     
     return(q)
}




