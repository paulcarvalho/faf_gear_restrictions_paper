### run_model
###
### Author: Paul Carvalho
###
### Description: Run fisheries model

run_model <- function(effort, gear.mgmt, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                      suit, ration, other, weight, sc_Linf, phi.min){
     ngears <- sum(gear.mgmt) # number of fishing gears
     N.ijte <- array(NA, dim = c(nsc, nspecies, t, length(effort)))  # empty shell to save data for scenario: abundance
     B.ijte <- array(NA, dim = c(nsc, nspecies, t, length(effort)))  # empty shell to save data for scenario: biomass
     cN.ijte <- array(NA, dim = c(nsc, nspecies, t, length(effort))) # empty shell to save data for scenario: catch (abundance)
     cB.ijte <- array(NA, dim = c(nsc, nspecies, t, length(effort))) # empty shell to save data for scenario: catch (biomass)
     S.ijte <- array(NA, dim = c(1, nspecies, t, length(effort)))    # empty shell to save data for scenario: spawning stock biomass
     R.ijte <- array(NA, dim = c(1, nspecies, t, length(effort)))    # empty shell to save data for scenario: recruitment
     
     for(e in 1:length(effort)){
        # determine fishing mortality given fishing effort
        fishing.effort <- (effort[e] / ngears) # divide effort among gears    
        fishing.effort <- fishing.effort * gear.mgmt
        F.mort <- F_mortality(nsc, nspecies, fishing.effort, gear.mgmt, q)
        
        # set up population
        N.ijt <- array(NA, dim = c(nsc, nspecies, t))
        N.ijt[,,1] <- N0
        S.ijt <- array(NA, dim = c(1, nspecies, t))
        R.ijt <- array(NA, dim = c(1, nspecies, t))
        c.ijt <- array(NA, dim = c(nsc, nspecies, t))
        c.ijt[,,1] <- 0
        # run model
        for(y in 2:t){
                # recruitment
                # calculate ssb and recruitment with Beverton-Holt
                S.ijt[,,y] <- calc_ssb(nspecies, nsc, L.lower, L.upper, N.ijt, y, Lmat)
                R.ijt[,,y] <- (S.ijt[,,y]) / (alpha + (beta*S.ijt[,,y]))
                N.tmp <- N.ijt[,,y-1]
                N.tmp[1,] <- N.tmp[1,] + R.ijt[,,y]
                
                # calculate predation mortality
                M2 <- calc_M2(N=N.tmp, suit, ration, nspecies, nsc, other, weight, sc_Linf, phi.min)
                
                # calculate catch
                c.tmp <- (F.mort / (F.mort + M1 + M2)) * ((1 - (exp(-(F.mort + M1 + M2)))) * N.tmp)
                
                # impose mortality
                N.tmp <- N.tmp * exp(-(F.mort + M1 + M2))
                
                # allow growth
                N.tmp <- calc_growth(N.tmp, phi, nsc, nspecies)
                # save abundance and catch
                N.ijt[,,y] <- N.tmp
                c.ijt[,,y] <- c.tmp
        }
        
        # calculate biomass and catch(biomass)
        B.ijt <- calc_bio(N.ijt, t, L.lower, L.upper, W.a, W.b, nsc, nspecies)
        cB.ijt <- calc_bio(c.ijt, t, L.lower, L.upper, W.a, W.b, nsc, nspecies)
        S.ijt <- calc_bio(S.ijt, t, L.lower, L.upper, W.a, W.b, 1, nspecies)
        
        # save data
        N.ijte[,,,e]  <- N.ijt
        B.ijte[,,,e]  <- B.ijt
        cN.ijte[,,,e] <- c.ijt
        cB.ijte[,,,e] <- cB.ijt
        S.ijte[,,,e]  <- S.ijt
        R.ijte[,,,e]  <- R.ijt
     }
     
     out <- list(N.ijte, B.ijte, cN.ijte, cB.ijte, S.ijte, R.ijte)
     return(out)
}