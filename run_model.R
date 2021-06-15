run_model <- function(effort, gear.mgmt, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min, N0){
     #' run_model
     #' 
     #' @description Run fisheries model with a single management scenario and set of parameters.
     #'
     #' @param effort Range of fishing effort
     #' @param gear.mgmt Management scenario
     #' @param nsc Number of size classes
     #' @param nspecies Number of functional groups, or species
     #' @param t Number of timesteps 
     #' @param Lmat Length at maturity
     #' @param M1 Natural mortality
     #' @param phi Proportion of fish from each functional group and size class that grow to the next size class
     #' @param L.lower Lower limit of each size class (cm)
     #' @param L.upper Upper limit of each size class (cm)
     #' @param W.a Length-weight conversion parameter
     #' @param W.b Length-weight conversion parameter
     #' @param q Product of catchability and selectivity
     #' @param alpha Beverton-Holt stock-recruitment parameter; 1/alpha is the max per capita production of recruits
     #' @param beta Beverton-Holt stock-recruitment parameter; recruitment approaches 1/beta as biomass increases
     #' @param suit Suitability of prey for predators
     #' @param ration Ration that must be consumed by predators to account for growth
     #' @param other Other food items that are not explicitly modeled
     #' @param weight Weight of each functional group and size class
     #' @param sc_Linf Size class at asymptotic length
     #' @param phi.min Minimum timestep
     #' @param N0 intial population size
     #'
     #' @return List with (1) abundance, (2) biomass, (3) catch in numbers, and (4) catch in biomass for each size class, functional 
     #' group, timestep, and effort level; and (5) spawning stock biomass and (6) recruitment for each functional group, timestep and
     #' effort level        

     ngears  <- sum(gear.mgmt) # number of fishing gears
     N.ijte  <- array(NA, dim = c(nsc, nspecies, t, length(effort)))  # empty shell to save data for scenario: abundance
     B.ijte  <- array(NA, dim = c(nsc, nspecies, t, length(effort)))  # empty shell to save data for scenario: biomass
     cN.ijte <- array(NA, dim = c(nsc, nspecies, t, length(effort))) # empty shell to save data for scenario: catch (abundance)
     cB.ijte <- array(NA, dim = c(nsc, nspecies, t, length(effort))) # empty shell to save data for scenario: catch (biomass)
     S.ijte  <- array(NA, dim = c(1, nspecies, t, length(effort)))    # empty shell to save data for scenario: spawning stock biomass
     R.ijte  <- array(NA, dim = c(1, nspecies, t, length(effort)))    # empty shell to save data for scenario: recruitment

     # Run fisheries model with range of effort levels in a parallel for loop
     n.cores <- parallel::detectCores() - 1
     cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
     registerDoParallel(cluster)

     results <- foreach(e = 1:length(effort), .combine = 'rbind', .multicombine = TRUE) %dopar% {
             source('model_functions.R')        
             # determine fishing mortality given fishing effort
             fishing.effort <- (effort[e] / ngears) # divide effort among gears
             fishing.effort <- fishing.effort * gear.mgmt
             F.mort         <- F_mortality(nsc, nspecies, fishing.effort, gear.mgmt, q)
             # set up population
             N.ijt      <- array(NA, dim = c(nsc, nspecies, t))
             N.ijt[,,1] <- N0
             S.ijt      <- array(NA, dim = c(1, nspecies, t))
             R.ijt      <- array(NA, dim = c(1, nspecies, t))
             c.ijt      <- array(NA, dim = c(nsc, nspecies, t))
             c.ijt[,,1] <- 0
             # run model
             for(y in 2:t){
                     # recruitment
                     # calculate ssb and recruitment with Beverton-Holt
                     S.ijt[,,y] <- calc_ssb(nspecies, nsc, L.lower, L.upper, N.ijt, y, Lmat)
                     R.ijt[,,y] <- (S.ijt[,,y]) / (alpha + (beta*S.ijt[,,y]))
                     N.tmp      <- N.ijt[,,y-1]
                     N.tmp[1,]  <- N.tmp[1,] + R.ijt[,,y]
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
             B.ijt  <- calc_bio(N.ijt, t, L.lower, L.upper, W.a, W.b, nsc, nspecies)
             cB.ijt <- calc_bio(c.ijt, t, L.lower, L.upper, W.a, W.b, nsc, nspecies)
             S.ijt  <- calc_bio(S.ijt, t, L.lower, L.upper, W.a, W.b, 1, nspecies)
             # save data
             list(N.ijt, B.ijt, c.ijt, cB.ijt, S.ijt, R.ijt)
     }
     parallel::stopCluster(cluster)
     
     for(i in 1:length(effort)){
             N.ijte[, , , i]  <- as.array(results[i, 1])[[1]]
             B.ijte[, , , i]  <- as.array(results[i, 2])[[1]]
             cN.ijte[, , , i] <- as.array(results[i, 3])[[1]]
             cB.ijte[, , , i] <- as.array(results[i, 4])[[1]]
             S.ijte[, , , i]  <- as.array(results[i, 5])[[1]]
             R.ijte[, , , i]  <- as.array(results[i, 6])[[1]]
     }
     
     out <- list(N.ijte, B.ijte, cN.ijte, cB.ijte, S.ijte, R.ijte)
     return(out)
}