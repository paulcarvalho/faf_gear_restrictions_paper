# Functions
# Author: Paul Carvalho

run_bsmodel <- function(nbs, landings, uvc, effort.bs, gear.mgmt, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min, N0){
  #' run_bsmodel
  #' @description Run the fisheries model and resample landings and UVC data via nonparametric bootstrapping to quantify uncertainty in catchability (function group/species) and selectivity (size).
  #'
  #' @param nbs number of bootstrapping iterations 
  #' @param landings fisheries-dependent data
  #' @param uvc fisheries-independent data
  #' @param effort.bs vector of effort levels to run 
  #' @param gear.mgmt see run_model() description
  #' @param nsc see run_model() description
  #' @param nspecies see run_model() description
  #' @param t see run_model() description
  #' @param Lmat see run_model() description
  #' @param M1 see run_model() description
  #' @param phi see run_model() description
  #' @param L.lower see run_model() description
  #' @param L.upper see run_model() description
  #' @param W.a see run_model() description
  #' @param W.b see run_model() description
  #' @param alpha see run_model() description
  #' @param beta see run_model() description
  #' @param suit see run_model() description
  #' @param ration see run_model() description
  #' @param other see run_model() description
  #' @param weight see run_model() description
  #' @param sc_Linf see run_model() description
  #' @param phi.min see run_model() description
  #' @param N0 see run_model() description
  #'
  #' @return List with (1) overall mean and confidence intervals for abundance, biomass, and catch and (2) mean and confidence intervals of abundance, biomass, and catch for each functional group/species.

  # nonparametric bootstrap resample landings and UVC data and rerun model
  num.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(my.cluster)
  
  results <- foreach(i = 1:nbs, .combine = 'rbind', .multicombine = TRUE) %dopar% {
    source('model_functions.R')
    source('run_model.R')
    library(dplyr)
    library(data.table)
    library(splitstackshape)
    library(plotrix)
    library(fitdistrplus)
    library(tidyr)
    library(foreach)
    library(doParallel)
    q.new <- NULL
    q.new <- bootstrap_qs(landings, uvc, resample = TRUE) # resample data via bootstrapping and generate bootstrapped catchability*selectivity parameters
    q.new <- array(as.numeric(unlist(q.new)), dim = c(nsc, nspecies, 3)) # turn list into array
    
    run.i <- run_model(effort = effort.bs, gear.mgmt = gear.mgmt, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q.new, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min, N0)
    # Total
    total.N  <- colSums(run.i[[1]][, , t, ], dims = 2)
    total.B  <- colSums(run.i[[2]][, , t, ], dims = 2)
    total.CN <- colSums(run.i[[3]][, , t, ], dims = 2)
    total.CB <- colSums(run.i[[4]][, , t, ], dims = 2)
    tmp1     <- data.frame(i = i, effort = effort.bs, Ntot = total.N, Btot = total.B, CNtot = total.CN, CBtot = total.CB)
    #tot.vals <- rbind(tot.vals, tmp1)
    # Functional group
    fg.N    <- as.vector(colSums(run.i[[1]][, , t, ]))
    fg.B    <- as.vector(colSums(run.i[[2]][, , t, ]))
    fg.CN   <- as.vector(colSums(run.i[[3]][, , t, ]))
    fg.CB   <- as.vector(colSums(run.i[[4]][, , t, ]))
    tmp2    <- data.frame(i = i, effort = rep(effort.bs, each = nspecies), fg = rep(c("Browser", "Detritivore", "Excavator/Scraper", "Grazer", "Macro-invertivore", "Micro-invertivore", "Pisci-invertivore", "Piscivore", "Planktivore"), 3), N = fg.N, B = fg.B, CN = fg.CN, CB = fg.CB)
    #fg.vals <- rbind(fg.vals, tmp2)
    
    list(tmp1, tmp2)
  }
  parallel::stopCluster(my.cluster)
  
  tmp3 <- NULL
  tmp4 <- NULL
  
  for(i in 1:nbs){
    tmp3 <- rbind(tmp3, results[i, 1][[1]]) 
    tmp4 <- rbind(tmp4, results[i, 2][[1]])
  }
  
  # calculate summary stats
  tmp3 <- tmp3 %>% dplyr::mutate(effort = as.factor(effort)) %>% dplyr::group_by(effort) %>% dplyr::summarise(N.mu = mean(Ntot), N.sd = sd(Ntot), B.mu = mean(Btot), B.sd = sd(Btot), CN.mu = mean(CNtot), CN.sd = sd(CNtot), CB.mu = mean(CBtot), CB.sd = sd(CBtot)) %>% dplyr::mutate(N.se = N.sd / sqrt(nbs), B.se = B.sd / sqrt(nbs), CN.se = CN.sd / sqrt(nbs), CB.se = CB.sd / sqrt(nbs)) %>% dplyr::mutate(N.lo = N.mu - qt(1 - (0.05 / 2), nbs - 1) * N.se, N.up = N.mu + qt(1 - (0.05 / 2), nbs - 1) * N.se, B.lo = B.mu - qt(1 - (0.05 / 2), nbs - 1) * B.se, B.up = B.mu + qt(1 - (0.05 / 2), nbs - 1) * B.se, CN.lo = CN.mu - qt(1 - (0.05 / 2), nbs - 1) * CN.se, CN.up = CN.mu + qt(1 - (0.05 / 2), nbs - 1) * CN.se, CB.lo = CB.mu - qt(1 - (0.05 / 2), nbs - 1) * CB.se, CB.up = CB.mu + qt(1 - (0.05 / 2), nbs - 1) * CB.se) %>% dplyr::select(effort, N.mu, N.up, N.lo, B.mu, B.up, B.lo, CN.mu, CN.up, CN.lo, CB.mu, CB.up, CB.lo)
  tmp4 <- tmp4 %>% dplyr::mutate(effort = as.factor(effort)) %>% dplyr::group_by(effort, fg) %>% dplyr::summarise(N.mu = mean(N), N.sd = sd(N), B.mu = mean(B), B.sd = sd(B), CN.mu = mean(CN), CN.sd = sd(CN), CB.mu = mean(CB), CB.sd = sd(CB)) %>% dplyr::mutate(N.se = N.sd / sqrt(nbs), B.se = B.sd / sqrt(nbs), CN.se = CN.sd / sqrt(nbs), CB.se = CB.sd / sqrt(nbs)) %>% dplyr::mutate(N.up = N.mu + qt(1 - (0.05 / 2), nbs - 1) * N.se, N.lo = N.mu - qt(1 - (0.05 / 2), nbs - 1) * N.se, B.up = B.mu + qt(1 - (0.05 / 2), nbs - 1) * B.se, B.lo = B.mu - qt(1 - (0.05 / 2), nbs - 1) * B.se, CN.up = CN.mu + qt(1 - (0.05 / 2), nbs - 1) * CN.se, CN.lo = CN.mu - qt(1 - (0.05 / 2), nbs - 1) * CN.se, CB.up = CB.mu + qt(1 - (0.05 / 2), nbs - 1) * CB.se, CB.lo = CB.mu - qt(1 - (0.05 / 2), nbs - 1) * CB.se)
  
  return(list(tmp3, tmp4))
}

calc_phi <- function(L.lower, L.upper, Linf, k, nsc){
  #' calc_phi
  #' 
  #' @description Calculates the proportion of fish from functional group i in size class j
  #' that grow to j+1 and re-scales times steps such that all fish in the fatest growing 
  #' functional group-size class comination reach the next size class. Code adapted from 
  #' Hall et al. (2006). 
  #'
  #' @param L.lower Lower limit of size class (included)
  #' @param L.upper Upper limit of size class (excluded)
  #' @param Linf Asymptotic length (cm)
  #' @param k Von Bertalanffy growth parameter
  #' @param nsc Number of size classes
  #'
  #' @return A list with (1) proportions of functional group i that grow from size class j to j+1 and (2) minimum time step
  
  # Upper and lower limits of size classes
  L.low <- repmat(L.lower, 1, length(Linf))
  L.up  <- repmat(L.upper, 1, length(Linf))
  # create matrix for VBGF parameters
  k.mat    <- repmat(t(as.matrix(k)), nsc, 1)
  Linf.mat <- repmat(t(as.matrix(Linf)), nsc, 1)
  # calculate time for fish to grow from lower to upper limit of size class
  tij           <- 1/k.mat * log((Linf.mat - L.low) / (Linf.mat - L.up)) # Hilborn and Walters 1992, p 428
  rownames(tij) <- as.character(L.lower)
  tij[which(is.infinite(tij))] <- 0 # replace Inf values with 0
  tij[which(is.na(tij))]       <- 0 # replace NaN with 0
  tij[which(tij < 0)]          <- 0 # replace negative values with 0
  
  # scale to the fastest growing species-size class combination
  phi_min <- min(tij[tij>0]) 
  for(i in 1:nspecies){
    for(j in 1:nsc){
      if(tij[j,i] > 0){
        tij[j,i] <- phi_min/tij[j,i]
      }
    }
  }
  return(list(tij, phi_min))
}

calc_growth <- function(N, phi, nsc, nspecies){
  #' calc_growth
  #'
  #' @description Calculate the number of fish that grow to next size class. Adopted from Hall et al. (2006).
  #'
  #' @param N Matrix with number of fish in size class j (row) and functional group i (col)
  #' @param phi Proportion of fish in size class j and functional group i that grow to the next size class 
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species, included in the model
  #'
  #' @return New abundance matrix (size class = row, functional group = col) after fish growth
  
  stay  <- N * (1-phi) 
  leave <- N * phi
  N.out <- stay + rbind(pracma::zeros(1, nspecies),leave[1:nsc-1,])
  return(N.out)
}

calc_ration <- function(k, Linf, nsc, nspecies, L.lower, L.upper, L.mid, W.a, W.b, phi.min, scale.Ge){
  #' calc_ration
  #' 
  #' @description Calculate the ration (ingestion rate) that must be consumed by species i in size class j 
  #' to account for modeled growth in a timestep (phi.min). Code adapted from Hall et al. (2006).
  #'
  #' @param k Von Bertalanffy growth parameter
  #' @param Linf Asymptotic length
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #' @param L.lower Lower limit of each size class (included)
  #' @param L.upper Upper limit of each size class (excluded)
  #' @param L.mid Midpoint of each size class
  #' @param W.a Length-weight conversion parameter
  #' @param W.b Length-weight conversion parameter
  #' @param phi.min Time step
  #' @param scale.Ge Scalar parameter for growth efficiency
  #'
  #' @return List with (1) ration that mught be consumed by each functional group-size class combination, (2) weight at midpoint of size class, 
  #' and (3) size class at asymptotic length (Linf)

  # Find the size class index at which each species/functional group reaches Linf
  sc_Linf <- rep(NA, nspecies)
  for(i in 1:nspecies){
    for(j in 1:nsc){
      if(Linf[i] > L.lower[j] & Linf[i] <= L.upper[j]){ sc_Linf[i] <- j }
        
    }
  }
  ration <- zeros(nrow = nsc, ncol = nspecies) # ration that must be consumed by each species/functional group-size class combination to account for growth.
  weight <- zeros(nrow = nsc, ncol = nspecies)
  g.eff  <- zeros(nrow = nsc, ncol = nspecies)  # Growth efficiency - proportion of food consumed that is converted to body mass
  
  for(i in 1:nspecies){
    for(j in 1:(sc_Linf[i]-1)){
      # Calculate the length after 1 timestep when the initial length is the mid point of a size class
      L2 <- L.mid[j] + (Linf[i] - L.mid[j]) * (1 - exp(-k[i] * phi.min))
      # Calculate growth increment from the weight at the mid-point
      weight[j,i] <- (W.a[i] * L.mid[j] ^ W.b[i])
      growth_inc  <- (W.a[i] * L2 ^ W.b[i]) - weight[j,i]
      g.eff[j,i]  <- 1 - (weight[j,i] / (W.a[i] * Linf[i] ^ W.b[i])) ^ 0.11 # equation 7 in Hall et al. (2006)
      g.eff[j,i]  <- g.eff[j,i] * 0.5 # 0.5 gives a reasonable scaling (Hall et al. 2006)
      g.eff[j,i]  <- g.eff[j,i] + (g.eff[j,i]*scale.Ge)
      ration[j,i] <- growth_inc * (1 / g.eff[j,i]) # ingestion rate that must be consumed by species/functional group i in size class j to account for modeled growth in a given timestep (tonnes/timestep)
    }
  }
  
  # repeat above for largest size class
  for(i in 1:nspecies){
    # Calculate the mid-point between the lower bound of the largest size class and the Linf
    fmid <- L.lower[sc_Linf[i]] + ((Linf[i] - L.lower[sc_Linf[i]]) / 2)
    # Calculate the length after 1 timestep when the initial length is the mid-point of a size class
    L2 <- fmid + (Linf[i] - fmid) * (1 - exp(-k[i] * phi.min))
    # Calculate growth increment from the weight at the mid-point
    weight[sc_Linf[i], i] <- W.a[i] * fmid ^ W.b[i]
    growth_inc            <- (W.a[i] * L2 ^ W.b[i]) - weight[sc_Linf[i], i]
    g.eff[sc_Linf[i], i]  <- 1 - (weight[sc_Linf[i],i] / (W.a[i] * Linf[i] ^ W.b[i])) ^ 0.11
    g.eff[sc_Linf[i], i]  <- g.eff[sc_Linf[i], i] * 0.5
    g.eff[sc_Linf[i], i]  <- g.eff[sc_Linf[i], i] + (g.eff[sc_Linf[i], i]*scale.Ge)
    ration[sc_Linf[i], i] <- growth_inc * (1 / g.eff[sc_Linf[i], i])
  }
  output <- list(ration, weight, sc_Linf)
  return(output) # Units:grams
}

calc_prefs <- function(L.mid, nsc, nspecies, mu, sigma, weight, sc_Linf){
  #' calc_prefs
  #' 
  #' @description Calculate the lognormal probability function for prey preferecnes, based on prey/predator
  #' size (i.e., weight) ratio. Code adapted from Hall et al. (2006).
  #'
  #' @param L.mid Midpoint of each size class (length in cm)
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #' @param mu Log(mean) for prey size size preference
  #' @param sigma Log(standard deviation) for prey size preference
  #' @param weight Weight
  #' @param sc_Linf Size class at asymptotic length
  #'
  #' @return Four-dimensional array with preferences for predator i in size class j and prey species l in size class k
  bs.ratio <- zeros(nsc, nspecies, nsc, nspecies)
  
  for(sp1 in 1:nspecies){ # predator species
    for(sp2 in 1:nspecies){ # prey species
      for(sc1 in 1:(sc_Linf[sp1])){ # only as far as the predators max size class
        for(sc2 in 1:(min(sc_Linf[sp1],sc_Linf[sp2]))){
          bs.ratio[sc1, sp1, sc2, sp2] <- (weight[sc2, sp2]) / (weight[sc1, sp1])
          bs.ratio[sc1, sp1, sc2, sp2] <- dlnorm(bs.ratio[sc1, sp1, sc2, sp2], mu, sigma)
        }
      }
    } 
  }
  
  return(bs.ratio)
}

calc_suit <- function(M2_prefs, tau, nsc, nspecies, sc_Linf){
  #' calc_suit
  #'
  #' @description Calculate the suitabilities for predator of size l, functional group k 
  #' for prey of size j, functional group i. Code adapted from Hall et al. (2006).
  #'
  #' @param M2_prefs Functional group and size class specific predator preferenes for prey
  #' @param tau Foodweb matrix
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #' @param sc_Linf Size class at asymptotic length for each functional group
  #'
  #' @return Four-dimensional matrix with suitabilities for predator and prey
  
  suit <- zeros(nsc, nspecies, nsc, nspecies)
  
  for(i in 1:nspecies){ # prey species
    for(j in 1:nsc){    # prey size class
      for(l in 1:nsc){  # predator size class
        suit[l,,j,i] <- M2_prefs[l,,j,i] * tau[,i] # equation 8.1 in Hall et al. (2006)
      }
    }
  }
  
  # standardize such that suitabilities sum to one
  # if the species is not a predator of anything, denom is zero
  denom <- matrix(NA, nsc, nspecies)
  for(sp1 in 1:nspecies){
    for(sc1 in 1:nsc){
      denom[sc1, sp1] <- sum(sum(suit[sc1,sp1,,]))
      if(denom[sc1, sp1] > 0){
        suit[sc1,sp1,,] <- suit[sc1,sp1,,] / denom[sc1, sp1]
      }
    }
  }
  return(suit)
}

calc_M2 <- function(N, suit, ration, nspecies, nsc, other, weight, sc_Linf, phi.min){
  #' calc_M2
  #' 
  #' @description Calculate the predation mortality (M2) for each functional group in each size class.
  #' Code adapted from Hall et al. (2006).
  #'
  #' @param N Abundance for size class j (rows) and functional group i (cols)
  #' @param suit Functional group-size class specific suitability of prey for predators
  #' @param ration Ingestion rate that accounts for fish growth
  #' @param nspecies Number of functional groups, or species
  #' @param nsc Number of size classes
  #' @param other Other food items that are not explicitly represented in the model
  #' @param weight Weight of each functional group and size class
  #' @param sc_Linf Size class at asymptotic length
  #' @param phi.min Time step scalar
  #'
  #' @return Matrix (size class = rows, functional group = cols) with predation mortality
  
  M2 <- matrix(0, nrow=nsc, ncol=nspecies)
  
  for(n in 1:nspecies){         # prey species
    for(m in 1:sc_Linf[n]){     # prey size class
      for(j in 1:nspecies){     # predator species
        for(i in 1:sc_Linf[j]){ # predator size class
          denom <- sum(suit[i,j,,] * weight * N) # denominator in equation 8 of Hall et al. (2006)
          denom <- denom + other
          if(denom > 0){
            M2[m,n] <- M2[m,n] + ration[i,j] * N[i,j] * (suit[i,j,m,n] / denom) # equation 8 of Hall et al. (2006)
          }
        }
      }
    }
  }
  
  output1 <- M2*phi.min # scale predation mortality to minimum time step (phi.min)
  return(output1)
}

nat_mortality <- function(L.lower, L.upper, nspecies, nsc, phi.min, Linf, k, which.L = 'mid'){
  #' nat_mortality
  #'
  #' @description Calculate natural mortality rates for each functional group and size class
  #' based on the equation from Gislason et al. (2010): ln(Mi,j) = 0.55 - 1.61 ln(Li,j) + 1.44 ln(Linf) + ln(ki).
  #'
  #' @param L.lower Lower limit of each size class 
  #' @param L.upper Upper limit of each size class
  #' @param nspecies Number of functional groups
  #' @param nsc Number of size classes
  #' @param phi.min Time step to scale rates
  #' @param Linf Asymptotic length
  #' @param k Von Bertalanffy growth parameter
  #' @param which.L Length used to estimate natural mortality ('lower', 'mid', 'upper', or 'mean'). Default = 'mid'
  #'
  #' @return Natural mortality rates for size class j (rows) and functional group i (cols) scaled to time step (phi.min)
  
  # create matrices for VBGF params
  Linf.mat <- repmat(t(as.matrix(Linf)), nsc, 1)
  k.mat    <- repmat(t(as.matrix(k)), nsc, 1)
  
  if(which.L == "mid"){
    # calculate midpoint of each size class
    Lmid     <- (L.lower + L.upper) / 2
    Lmid.mat <- repmat(as.matrix(Lmid), 1, nspecies)
    # calculate natural mortality and scale to time step of fastest growing species-size class combination
    lnM    <- 0.55 - (1.61*log(Lmid.mat)) + (1.44*log(Linf.mat)) + (log(k.mat))
    M.year <- exp(lnM)
    M      <- M.year*phi.min # scale to timestep
    # Plot natural mortality curves
    cbp2   <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    M.df   <- as.data.frame(M)
    M.plot <- M.df %>%
      dplyr::mutate(size = Lmid) %>%
      dplyr::select(size, "Browsers" = V1, "Detritivores" = V2, "Excavators/scrapers" = V3,
                    "Grazers" = V4, "Macro-invertivores" = V5, "Micro-invertivores" = V6,
                    "Pisci-invertivores" = V7, "Piscivores" = V8, "Planktivores" = V9) %>%
      gather(., "fg", "M", -size) %>%
      filter(fg != "Browsers" | size <= 37.5) %>%
      filter(fg != "Detritivores" | size <= 17.5) %>%
      filter(fg != "Excavators/scrapers" | size <= 27.5) %>%
      filter(fg != "Grazers" | size <= 37.5) %>%
      filter(fg != "Macro-invertivores" | size <= 37.5) %>%
      filter(fg != "Micro-invertivores" | size <= 27.5) %>%
      filter(fg != "Planktivores" | size <= 37.5) %>%
      ggplot() +
        geom_line(aes(x = size, y = M, color = fg), lwd = 1.5, alpha = 0.5) +
        theme_classic() +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        scale_color_manual(values = cbp2) +
        theme(legend.title = element_blank(),
              legend.position = c(0.8,0.7)) +
        labs(x = "Size (cm)", y = "Natural mortality (M/timestep)")
    print(M.plot)
    
  } else if (which.L == "upper") {
    # create matrix
    Lup.mat <- repmat(as.matrix(L.upper), 1, nspecies)
    # calculate natural mortality and scale to time step of fastest growing species-size class combination
    lnM    <- 0.55 - (1.61*log(Lup.mat)) + (1.44*log(Linf.mat)) + (log(k.mat))
    M.year <- exp(lnM)
    M      <- M.year*phi.min # scale to timestep
    
  } else if(which.L == "lower") {
    # create matrix
    Llo.mat <- repmat(as.matrix(L.lower), 1, nspecies)
    # calculate natural mortality and scale to time step of fastest growing species-size class combination
    lnM    <- 0.55 - (1.61*log(Llo.mat)) + (1.44*log(Linf.mat)) + (log(k.mat))
    M.year <- exp(lnM)
    M      <- M.year*phi.min # scale to timestep
    
  } else if(which.L == "mean") {
    L.add <- 0:((L.upper[1] - L.lower[1]) - 1)
    M.tmp <- array(NA, dim = c(nsc, nspecies, length(L.add)))
    for(i in 1:length(L.add)){
      L          <- L.lower + L.add[i]
      L          <- repmat(as.matrix(L), 1, nspecies)
      lnM        <- 0.55 - (1.61*log(L)) + (1.44*log(Linf.mat)) + (log(k.mat))
      M.year     <- exp(lnM)
      M.tmp1     <- M.year*phi.min # scale to timestep
      M.tmp[,,i] <- M.tmp1
    }
    M <- apply(M.tmp, c(1,2), mean)
  }
  return(M)
}

calc_N0 <- function(uvc){
  #' calc_N0
  #'
  #' @param uvc Empirical underwater census (fisheries-independent) data
  #'
  #' @return Matrix with initial starting abundance for each functional group (or species) and size class

  library(DescTools)
  
  # Calculate mean, standard deviation, and standard error for abundance per hectare
  fg.abundance <- uvc %>% dplyr::group_by(site_name, transect, fg) %>% dplyr::summarise(abundance = sum(abundance)) %>% dplyr::group_by(fg) %>% dplyr::summarise(mean_ab = mean(abundance), sd_ab = sd(abundance), se_ab = std_err(abundance)) %>% dplyr::mutate(mean_ab = mean_ab *40, sd_ab = sd_ab * 40, se_ab = se_ab * 40)
  
  # Size spectra functions to calculate N0
  inverse_method <- function(n, xmin, xmax, b){
    u <- runif(round(n))
    x <- ((u*xmax^(b+1)) + ((1-u)*xmin^(b+1)))^(1/(b+1))
    return(x)
  }
  bin_sizes <- function(n, xmax){
    nmax <- RoundTo(max(n),(maxsize-minsize)/nsc,floor)
    x <- cut(n, breaks=c(seq(minsize,nmax+(maxsize-minsize)/nsc,(maxsize-minsize)/nsc)), 
             labels=as.character(c(seq(minsize,nmax,(maxsize-minsize)/nsc))), include.lowest = TRUE, right = FALSE)
    if(nmax < xmax){
      levels(x) = as.character(c(seq(minsize,xmax,(maxsize-minsize)/nsc)))
    }
    mat <- as.matrix(summary(x))
    return(mat)
  }
  
  # Length spectra values. Values fall within estimates derived in Robinson et al. (2017).
  # Herbivores and detritivores were set shallower than other functional groups due to
  # metabolic constraints on abundance-size relationship (Trebilco et al. 2013,
  # Robinson and Baum 2019).
  # b = -2.3 # herbivores and detritivores
  # b = -2.7 # piscivores and invertivores and planktivores
  fg <- c("Browser", "Detritivore", "Excavator/scraper", "Grazer", "Macro-invertivore", "Micro-invertivore", "Pisci-invertivore", "Piscivore", "Planktivore")
  
  N0 <- matrix(NA, nrow = nsc, ncol = length(Linf), dimnames = list(as.character(L.lower), NULL))
  for(i in 1:length(Linf)){
    if(fg[i] == "Browser" | fg[i] == "Excavator/scraper" |
       fg[i] == "Grazer" | fg[i] == "Detritivore") b = -2.3 else b = -2.7
       abundance <- fg.abundance %>% filter(fg == fg[i]) %>% dplyr::select(mean_ab)
       n <- inverse_method(abundance$mean_ab, Lmin[i], Lmax[i], b)
       Nij <- bin_sizes(n, Lmax[i])
       if(length(Nij) > nsc) Nij = Nij[1:nsc,]
       tmp.mat <- matrix(0, nrow = nsc, ncol = 1)
       tmp.mat[1:length(Nij),1] <- Nij
       N0[,i] <- tmp.mat
  }
  rownames(N0) <- as.character(L.lower)
  return(N0)
}

bootstrap_qs <- function(landings, uvc, resample = NULL){
  #' bootstrap_qs
  #'
  #' @description Resample (via bootstrapping) landings and uvc data to quantify uncertainty in catchability and selectivity.
  #'
  #' @param landings Empirical landings (fisheries-dependent) data
  #' @param uvc Empirical underwater census (fisheries-independent) data  
  #' @param resample TRUE to bootstrap; FALSE to use original dataset
  #'
  #' @return List with products of catchability and selectivity for (1) hook, (2) net, and (3) spear.

  # Filter out size that are not included in the model and not exploited in the fishery
  landings.tmp <- landings %>% filter(size_cm >= 10 & size_cm < 65)
  uvc.tmp      <- uvc %>% filter(size_cm >= 10 & size_cm < 65) %>% dplyr::mutate(site_tran = paste(site_name, transect))
  uni.fg       <- unique(landings.tmp$fg) # Get unique functional groups
  
  if(resample == FALSE){
    uvc.tmp <- uvc.tmp %>% dplyr::select(fg, size_5cm_bin, site = site_name, site_tran, abundance)
    uvc.tmp <- setDT(expandRows(uvc.tmp, 'abundance'))
  } else if(resample == TRUE){
    # Resample landings data via bootstrapping
    trip.sample <- landings.tmp %>% # Save trip IDs and gear and number of fish caught
                   dplyr::group_by(trip_id, fishing_ground, gear_cat1) %>% 
                   dplyr::summarize(size = sum(abundance)) 
    landings.bs <- NULL # empty dataframe to fill during resampling of landings data
    for(i in 1:length(trip.sample$trip_id)){             # iterate through trips and resample from all catch from each fishing ground
      trip.i           <- trip.sample$trip_id[i]         # index for a single trip
      gear.i           <- trip.sample$gear_cat1[i]       # get gear for trip.i
      fish.grnd.i      <- trip.sample$fishing_ground[i]  # get fishing ground for trip.i
      sample.i         <- landings.tmp %>%         # save fg and sizes caught for trip.i
                          # dplyr::filter(trip_id == trip.i) %>%
                          dplyr::filter(fishing_ground == fish.grnd.i & gear_cat1 == gear.i) %>%
                          dplyr::select(fg, size_cm, bin_5cm) %>% 
                          filter(size_cm >= 10) 
      index.bs         <- sample(x = 1:length(sample.i$fg), size = trip.sample$size[i], replace = TRUE) # sample from trip.i with replacement
      tmp.bs           <- sample.i[index.bs,] # get functional groups and size classes for new sample
      tmp.bs$trip_id   <- trip.i
      tmp.bs$gear_cat1 <- gear.i
      landings.bs      <- rbind(landings.bs, tmp.bs)
    }
    landings.tmp <- landings.bs
    
    # Resample UVC data via bootstrapping
    site.sample <- uvc.tmp %>% dplyr::mutate(site_tran = paste(site_name, transect)) %>% dplyr::select(site_name, site_tran, fg, size_5cm_bin, abundance) %>% dplyr::group_by(site_name, site_tran) %>% dplyr::summarize(size = sum(abundance)) # Save site and transects and number of fish observed
    uvc.bs      <- NULL # empty dataframe to fill during resampling of landings data
    for(i in 1:length(site.sample$site_tran)){     # Iterate through site_trans and resample from fish observed at each site
      site.i           <- site.sample$site_name[i] # save site name for site_tran[i]
      sample.i         <- uvc.tmp %>%              # save all fish observed at site.i
                          dplyr::filter(site_name == site.i) %>% 
                          dplyr::select(fg, size_5cm_bin, abundance) 
      sample.i         <- setDT(expandRows(sample.i, 'abundance')) # expand table based on abundance
      index.bs         <- sample(x = 1:length(sample.i$fg), size = site.sample$size[i], replace = TRUE) # sample indexes from all fish observed at site.i
      tmp.bs           <- sample.i[index.bs, ] # save the functional groups and size classes that were sampled
      tmp.bs$site      <- site.i # save site name
      tmp.bs$site_tran <- site.sample$site_tran[i] # save the site name and transect
      uvc.bs           <- rbind(uvc.bs, tmp.bs) # add to full bootstrapped sample
    }
    uvc.tmp <- uvc.bs
  }
  
  # Calculate the ratio of catch per trip for each functional group (i), and functional group and size class (ij) - bootstrapped landings data
  c_i         <- landings.tmp %>% 
                 dplyr::group_by(trip_id, gear_cat1, fg) %>% 
                 dplyr::mutate(abundance = 1) %>% 
                 dplyr::summarise(abundance = sum(abundance)) %>% 
                 dplyr::group_by(trip_id, gear_cat1) %>% 
                 dplyr::mutate(ratio = abundance/sum(abundance)) %>% 
                 dplyr::group_by(fg, gear_cat1) %>% 
                 dplyr::summarise(c = mean(ratio), c_stderr = std.error(ratio))
  names(c_i)  <- c("fg", "gear", "c", "c_stderr")
  c_ij        <- landings.tmp %>% 
                 dplyr::group_by(trip_id, gear_cat1, fg, bin_5cm) %>% 
                 dplyr::mutate(abundance = 1) %>% 
                 dplyr::summarise(abundance = sum(abundance)) %>% 
                 dplyr::group_by(trip_id, gear_cat1) %>% 
                 dplyr::mutate(ratio = abundance/sum(abundance)) %>%
                 dplyr::group_by(fg, bin_5cm, gear_cat1) %>% 
                 dplyr::summarise(c = mean(ratio), c_stderr = std.error(ratio))
  names(c_ij) <- c("fg", "sc", "gear", "c", "c_stderr")
  
  # Calculate the ratio of numbers per transect per functional group (i) and functional group and size class (ij) - bootstrapped UVC data  
  n_i         <- uvc.tmp %>% 
                 dplyr::group_by(site_tran, fg) %>% 
                 dplyr::mutate(abundance = 1) %>% 
                 dplyr::summarise(abundance = sum(abundance)) %>% 
                 dplyr::group_by(site_tran) %>% 
                 dplyr::mutate(ratio = abundance/sum(abundance)) %>% 
                 dplyr::group_by(fg) %>% 
                 dplyr::summarise(n = mean(ratio), n_stderr = std.error(ratio))
  names(n_i)  <- c("fg", "n", "n_stderr")
  n_ij        <- uvc.tmp %>% 
                 dplyr::group_by(site_tran, fg, size_5cm_bin) %>% 
                 dplyr::mutate(abundance = 1) %>% 
                 dplyr::summarise(abundance = sum(abundance)) %>%
                 dplyr::group_by(site_tran) %>% 
                 dplyr::mutate(ratio = abundance/sum(abundance)) %>% 
                 dplyr::group_by(fg, size_5cm_bin) %>% 
                 dplyr::summarise(n = mean(ratio), n_stderr = std.error(ratio))
  names(n_ij) <- c("fg", "sc", "n", "n_stderr")  
  
  cn_ij   <- c_ij %>% full_join(., n_ij, by = c("fg", "sc")) # Merge landings and UVC dataframes
  q_i     <- c_i %>% full_join(., n_i, by = c("fg")) %>% dplyr::mutate(q = c/n) # Calculate catchability; fisheries-dependent/-independent
  L.lower <- seq(from = 0, to = 65 - (65-5)/12, length.out = 13)   # low limit included in the size class
  L.upper <- seq(from = 0 + ((65-5)/12), to = 65, length.out = 13) # upper limit excluded
  Lmid    <- ((L.lower + L.upper) / 2)[2:12] # Only include exploited size classes or modeled size classes
  
  hook.dist  <- landings.tmp %>% filter(gear_cat1 == "line");  hf.lnorm <- fitdist(hook.dist$size_cm, "lnorm") # fit lognormal distribution to hook-and-line landings
  net.dist   <- landings.tmp %>% filter(gear_cat1 == "net");   nf.lnorm <- fitdist(net.dist$size_cm, "lnorm") # fit lognormal distribution to net landings data
  spear.dist <- landings.tmp %>% filter(gear_cat1 == "spear"); sf.lnorm <- fitdist(spear.dist$size_cm, "lnorm") # fit lognormal distribution to spear landings data
  
  # Optimization function
  q_optim <- function(params, x, mu, sd, cn_ij, gear_type, fg, sc){
    q <- params[1:9] # functional group catchabilities, calibrated through optimization
    # model
    sj    <- dlnorm(x, meanlog = mu, sdlog = sd)
    sj    <- as.matrix(sj)
    model <- sj %*% q
    # expected
    cn_ij_x    <- cn_ij %>% dplyr::mutate(fg = as.factor(fg)) %>% dplyr::filter(gear == gear_type)
    cn_ij_x$cn <- cn_ij_x$c / cn_ij_x$n
    cn_ij_x    <- cn_ij_x[,c(1,2,8)]
    ## need to check for missing fg and sc combinations
    tmp1    <- paste(cn_ij_x$fg, cn_ij_x$sc, sep=' ')
    uni.tmp <- unique(tmp1)
    tmp2    <- paste(rep(uni.fg, each=length(sc)), rep(sc, times=length(uni.fg)), sep = ' ')
    tmp3    <- setdiff(tmp2, tmp1)
    tmp4    <- data.frame(fg = sapply(strsplit(tmp3, split = ' '), '[', 1), sc = sapply(strsplit(tmp3, split = ' '), '[', 2), cn = NA)
    cn_ij_x <- rbind(cn_ij_x, tmp4)
    cn_ij_x_wide <- spread(cn_ij_x, fg, cn, drop = FALSE)
    cn_ij_x_mat  <- as.matrix(cn_ij_x_wide[,-1])
    expected     <- cn_ij_x_mat
    # calculate sum of squared errors
    sse <- sum(((model) - (expected))^2, na.rm = TRUE)
    return(sse)
  }
  
  # Hook-and-line catchability (q) optimization
  par      <- subset(q_i, gear == "line")[,c(1,7)] # initial catchabilities for line fishing
  missing  <- setdiff(uni.fg, par$fg) # check if functional groups are missing
  if (length(missing) > 0) {par <- rbind(data.frame(fg = missing, q = 0), (as.data.frame(par)))} # fill in missing functional group(s) with q=0
  hook_par <- t(as.matrix(par$q))  # save catchabilities for optimization
  mu.hook  <- hf.lnorm$estimate[1] # log(mean) for size selectivity
  sd.hook  <- hf.lnorm$estimate[2] # log(standard deviation) for size selectivity
  x_hook   <- Lmid # Set size classes to model
  hook_fit <- optim(par = hook_par, fn = q_optim, lower = c(0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), upper = c(1, 500, 500, 500, 500, 500, 500, 500, 500), method = "L-BFGS-B", control = list(maxit = 5000), x = x_hook, mu = mu.hook, sd = sd.hook, cn_ij = cn_ij, gear_type = "line", fg = unique(landings.tmp$fg), sc = unique(landings.tmp$bin_5cm))
  F_hook             <- (as.matrix(dlnorm(x_hook,mu.hook,sd.hook))) %*% (as.matrix(hook_fit$par)) # Calculate q*s with optimized parameters
  F_hook[F_hook < 0] <- 0 # Set negative values, if any, to zero
  F_hook             <- F_hook / sum(F_hook) # Calculate proportional q*s values
  
  # Net catchability (q) optimization
  par     <- subset(q_i, gear == "net")[,c(1,7)] # initial catchabilities for net fishing
  missing <- setdiff(uni.fg, par$fg) # check if functional groups are missing
  if (length(missing) > 0) {par <- rbind(data.frame(fg = missing, q = 0), (as.data.frame(par)))} # fill in missing functional group(s) with q=0
  net_par <- t(as.matrix(par$q)) # save catchabilities for optimization
  mu.net  <- nf.lnorm$estimate[1] # log(mean) for size selectivity
  sd.net  <- nf.lnorm$estimate[2] # log(standard deviation) for size selectivity
  x_net   <- Lmid # Set size classes to model
  net_fit <- optim(par = net_par, fn = q_optim, lower = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), upper = c(500, 500, 500, 500, 500, 500, 500, 500, 500), method = "L-BFGS-B", control = list(maxit = 10000), x = x_net, mu = mu.net, sd = sd.net, cn_ij = cn_ij, gear_type = "net", fg = unique(landings.tmp$fg), sc = unique(landings.tmp$bin_5cm))
  F_net            <- (as.matrix(dlnorm(x_net, mu.net, sd.net))) %*% (as.matrix(net_fit$par)) # Calculate q*s with optimized parameters
  F_net[F_net < 0] <- 0 # Set negative values, in any, to zero
  F_net            <- F_net / sum(F_net) # Caclulate proportional q*s values
  
  # Spear catchability (q) optimization
  par       <- subset(q_i, gear == "spear")[,c(1,7)] # initial catchabilities for spear fishing
  missing   <- setdiff(uni.fg, par$fg) # check if functional groups are missing
  if (length(missing) > 0) {par <- rbind(data.frame(fg = missing, q = 0), (as.data.frame(par)))} # fill in missing functional group(s) with q=0
  spear_par <- t(as.matrix(par$q)) # save catchabilities for optimization
  mu.spear  <- sf.lnorm$estimate[1] # log(mean) for size selectivity
  sd.spear  <- sf.lnorm$estimate[2] # log(standard deviation) for size selectivity
  x_spear   <- Lmid # Set size classes to model  
  spear_fit <- optim(par = spear_par, fn = q_optim, lower = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), upper = c(500, 500, 500, 500, 500, 500, 500, 500, 500), method = "L-BFGS-B", control = list(maxit = 10000), x = x_spear, mu = mu.spear, sd = sd.spear, cn_ij = cn_ij, gear_type = "spear", fg = unique(landings.tmp$fg), sc = unique(landings.tmp$bin_5cm))
  F_spear              <- (as.matrix(dlnorm(x_spear, mu.spear, sd.spear))) %*% (as.matrix(spear_fit$par)) # Calculate q*s with optimized parameters
  F_spear[F_spear < 0] <- 0 # Set negative values, if any, to zero
  F_spear              <- F_spear / sum(F_spear) # Calculate proportional q*s values
  
  # PLOTS
  # hook_heat  <- heatmap(t(F_hook), Colv=NA, Rowv=NA, scale="none", labCol=levels(cn_ij$sc), labRow=levels(as.factor(cn_ij$fg)), margins=c(9.5,9.5))
  # net_heat   <- heatmap(t(F_net), Colv=NA, Rowv=NA, scale="none", labCol=levels(cn_ij$sc), labRow=levels(as.factor(cn_ij$fg)), margins=c(9.5,9.5))
  # spear_heat <- heatmap(t(F_spear), Colv=NA, Rowv=NA, scale="none", labCol=levels(cn_ij$sc), labRow=levels(as.factor(cn_ij$fg)), margins=c(9.5,9.5))
  
  return(list(F_hook, F_net, F_spear))  
}

F_mortality <- function(nsc, nspecies, fishing.effort, gear.mgmt, q){
  #' F_mortality
  #'
  #' @description Calculate fishing mortality for each gear.
  #'
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #' @param fishing.effort Total fishing effort (arbitrary value)
  #' @param gear.mgmt Management strategy to simulate
  #' @param q Product of catchability and selectivity for each gear type
  #'
  #' @return Fishing mortalities
  
  F.tmp <- array(NA, dim = c(nsc, nspecies, length(gear.mgmt))) # empty shell for fishing mortality of each fishing gear
  F.out <- array(0, dim = c(nsc, nspecies))
  # set up gear management scenario
  q.tmp <- array(NA, dim = c(nsc, nspecies, length(gear.mgmt)))
  for(i in 1:length(gear.mgmt)){
    if(gear.mgmt[i] == 1){
      q.tmp[,,i] <- q[,,i]/sum(q[,,i]) # assign catchability/selectivity and scale such that qs for all sums to one
    } else {
      q.tmp[,,i] <- 0
    }
  }
  # calculate fishing mortality
  for(i in 1:length(gear.mgmt)){
    F.tmp[,,i] <- fishing.effort[i] * q.tmp[,,i]
    F.out      <- F.out + F.tmp[,,i]
  }
  return(F.out) # I don't think necessary to scale with time step
}

calc_catch <- function(N, F.mort){
  #' calc_catch
  #'
  #' @description Calculate the catch in numbers.
  #'
  #' @param N Abundance of each functional group, or species, in each siez class
  #' @param F.mort Fishing mortliaty
  #'
  #' @return Total catch
  
  tot.catch <- N * (1 - exp(-F.mort)) # calculate F mortality
  # Mathimatically the following won't occur, but this was added for old calculation of fishing mortliaty and kept
  if(sum(c(N - tot.catch) < 0) > 0) stop("Error: catch for one or more species-size class combinations exceeds limit")
  return(tot.catch)
}

calc_bio <- function(N, t, L.lower, L.upper, W.a, W.b, nsc, nspecies){
  #' calc_bio
  #'
  #' @description Convert length to biomass using the midpoint of each size class.
  #'
  #' @param N Abundance of each functional group and size class 
  #' @param t Number of timesteps
  #' @param L.lower Lower limit of each size class
  #' @param L.upper Upper limit of each size class
  #' @param W.a Length-weight conversion parameter
  #' @param W.b Length-weight conversion parameter
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #'
  #' @return Biomass of each functional group and size class
  
  # Save N as an array if it's a matrix
  if(t == 1){ N <- array(N, dim=c(nsc,nspecies,t)) }
  # calculate midpoint of each size class
  Lmid       <- (L.lower + L.upper) / 2
  Lmid.mat   <- pracma::repmat(as.matrix(Lmid), 1, nspecies)     
  Lmid.array <- array(Lmid.mat, dim=c(nsc,nspecies,t))
  # create matrices for params
  wa.mat   <- pracma::repmat(t(as.matrix(W.a)), nsc, 1)
  wa.array <- array(wa.mat, dim=c(nsc,nspecies,t))
  wb.mat   <- pracma::repmat(t(as.matrix(W.b)), nsc, 1)
  wb.array <- array(wb.mat, dim=c(nsc,nspecies,t))
  # calculate biomass (kg)
  N.out <- N * ((wa.array*Lmid.array^wb.array) / 1000)
  return(N.out)
}

calc_propFG <- function(N_bio){
  #' calc_propFG
  #'
  #' @description Calculate the proportion of biomass in each functional group, or species
  #'
  #' @param N_bio Biomass in each functional group and size class
  #'
  #' @return Proportion of total biomass in each functional group
  sum_bio  <- colSums(N_bio)
  prop_bio <- sum_bio / sum(sum_bio)
  return(prop_bio)
}

calc_fgsize_output <- function(X, effort, total.effort, nsc, nspecies){
  #' calc_fgsize_output
  #' 
  #' @description Calculate the biomass or catch and size distribution for each functional group and effort combination.
  #'
  #' @param X Biomass or catch
  #' @param effort Range of fishing effort
  #' @param total.effort Total fishing effort
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional groups, or species
  #'
  #' @return Proportion of fish in each size class for each functional group
  
  # empty dataframes for saving data
  out.df  <- data.frame(fg = as.character(), E = as.numeric(), V = as.numeric(), V.rel = as.numeric(), sizedist = as.numeric())
  out.df2 <- data.frame(E = as.numeric(), V = as.numeric())
  func.groups <- c("Browser","Detritivore","Excavator/scraper","Grazer","Macro-invertivore","Micro-invertivore","Pisci-invertivore","Piscivore","Planktivore")
  for(i in 1:length(effort)){
    effort.tmp <- rep((effort[i] / total.effort), nspecies)
    X.tmp      <- X[,,i]
    # calculate proportions in each size class
    X.sum  <- repmat(t(as.matrix(colSums(X.tmp))), nsc, 1)
    X.prop <- X.tmp / X.sum
    # get only proportion in the smallest two size classes, that is
    # the size class of recruits and smallest size class exposed to fishing
    X.prop <- colSums(X.prop[1:2,])
    # calculate total biomass
    V.tmp  <- colSums(X.tmp)
    V2.tmp <- sum(V.tmp)
    # calculate relative biomass for each functional group
    V.rel.tmp <- V.tmp / sum(V.tmp)
    # save all data
    out.tmp  <- data.frame(fg=func.groups, E = effort.tmp, V = V.tmp, V.rel = V.rel.tmp, sizedist = X.prop)
    out.tmp2 <- data.frame(E = effort.tmp[1], V = V2.tmp)
    # combine with other data
    out.df  <- rbind(out.df, out.tmp)
    out.df2 <- rbind(out.df2, out.tmp2)
  }
  return(list(out.df, out.df2))
}

cal_recruitment <- function(par.in, N0, t, nsc, nspecies, M1, phi, L.lower, L.upper, W.a, W.b, Bi.F0, suit, ration, weight, sc_Linf, phi.min){
  #' calc_recruitment
  #' 
  #' @description Calibrate recruitment such that biomass of each functional group, or species, reaches biomass/ha
  #' in remote areas of Indonesia (Campbell et al. 2020). The main script calls this function within optim to find
  #' recruitment that minimizes sum of squared residuals.
  #'
  #' @param par.in (1) recruitment rate and (2) other food items
  #' @param N0 Initial starting population
  #' @param t Number of timesteps
  #' @param nsc Number of size classes
  #' @param nspecies Number of functional gropus, or species
  #' @param M1 Natural mortality
  #' @param phi Proportion of fish in size class j and functional group i that grow to the next size class
  #' @param L.lower Lower limit of each size class (cm)
  #' @param L.upper Upper limit of each size class (cm)
  #' @param W.a Length-weight conversion parameter
  #' @param W.b Length-weight conversion parameter
  #' @param Bi.F0 Expected pristine biomass (i.e., at zero fishing effort)
  #' @param suit Suitability of prey for predators
  #' @param ration Ration that must be consumed to account for growth
  #' @param weight Weight of each functional group and size class
  #' @param sc_Linf Size class at asymptotic length
  #' @param phi.min Minimum timestep
  #'
  #' @return Residuals of expected pristine biomass and model biomass at the last timestep
  
  r.i   <- par.in[1:9]
  other <- par.in[10]
  N.ijt <- array(NA, dim = c(nsc, nspecies, t)) # empty shell to save data
  N.ijt[,,1] <- N0 # set initial abundance
  # run model
  for(i in 2:t){
    # recruitment
    N.tmp     <- N.ijt[,,i-1]
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

calc_ssb <- function(nspecies, nsc, L.lower, L.upper, N.ijt, y, Lmat){
  #' calc_ssb
  #' 
  #' @description Calculate the spawning stock biomass for each functional group, or species.
  #'
  #' @param nspecies Number of functional groups, or species
  #' @param nsc Number of size classes
  #' @param L.lower Lower limit of size class (cm)
  #' @param L.upper Upper limit of size class (cm)
  #' @param N.ijt Abundance for each functional group in each size class for all timesteps
  #' @param y Year index
  #' @param Lmat Length at maturity 
  #'
  #' @return Spawning stock biomass of each functional group
  
  # calc spawning stock biomass
  ssb <- vector("numeric", length = nspecies)
  for(x in 1:nspecies){
    stock <- NULL
    for(z in 1:nsc){
      tmp_stock <- data.frame(size = seq(L.lower[z]+1,L.upper[z],by=1), stock = N.ijt[z,x,y-1]/5)
      stock     <- rbind(stock, tmp_stock)
    }
    ssb[x] <- sum(subset(stock, size >= Lmat[x])$stock)
  }   
  return(ssb)
}

calc_sensitivity_indices <- function(N.ijte, B.ijte, cN.ijte, cB.ijte, L.mid, Lmat, tau, base.effort){
  #' calc_sensitivity_indices
  #' 
  #' @description Calculate indicators for sensitiviy to model parameters
  #'
  #' @param N.ijte Abundance
  #' @param B.ijte Biomass
  #' @param cN.ijte Catch in numbers
  #' @param cB.ijte Catch in biomass
  #' @param L.mid Midpoint of each size class
  #' @param Lmat Length at maturity
  #' @param tau Food web matrix
  #' @param base.effort Original fishing effort 
  #'
  #' @return Dataframe with effort, biomass, catch in biomass, biomass size distribution, and catch size distribution
  
  # Calculate total biomass at equilibrium for all efforts
  B.equil <- colSums(colSums(B.ijte[,,t,]))
  # Calculate total catch at equilibrium for all efforts
  cB.equil <- colSums(colSums(cB.ijte[,,t,]))
  # Calculate mean length of stock and catch
  BmeanL <- vector(mode="numeric", length=length(effort))
  CmeanL <- vector(mode="numeric", length=length(effort))
  for(e in 1:length(base.effort)){
    totalN    <- sum(N.ijte[,,t,e])
    totalL    <- sum(rowSums(N.ijte[,,t,e]) * L.mid)
    BmeanL[e] <- totalL/totalN
    totalcN   <- sum(cN.ijte[,,t,e])
    totalcL   <- sum(rowSums(cN.ijte[,,t,e]) * L.mid)
    CmeanL[e] <- totalcL/totalcN
  }
  # Create dataframe
  si.df <- data.frame(effort = base.effort, B = B.equil, cB = cB.equil, BmeanL = BmeanL, CmeanL = CmeanL)
  # Find values at 0.5B0
  si.df <- si.df %>% dplyr::mutate(sq.diff = (B - (max(B)/2))^2)
  out   <- si.df[which(si.df$sq.diff == min(si.df$sq.diff)),]
  return(out)
}

calc_summary_indices <- function(N.ijte, B.ijte, cN.ijte, cB.ijte, t, L.mid, Lmat, nspecies){
  #' calc_summary_indices
  #' 
  #' @description Calculate summary for model outputs
  #'
  #' @param N.ijte Abundance
  #' @param B.ijte Biomass
  #' @param cN.ijte Catch in numbers
  #' @param cB.ijte Catch in biomass
  #' @param t Number of timesteps
  #' @param L.mid Midpoint of size classes
  #' @param Lmat Length at maturity
  #' @param nspecies Number of functional groups, or species
  #'
  #' @return List with (1) total biomass, (2) total catch biomass, (3) biomass of each functional group, (4) catch biomass of each functional group,
  #' (5) mean weight, (6) mean length, (7) mean biomass at length at maturity, (8) mean weight of catch, (9) and mean length of catch
  
  B.equil    <- colSums(colSums(B.ijte[,,t,])) # Calculate total biomass at equilibrium for all efforts
  cB.equil   <- colSums(colSums(cB.ijte[,,t,])) # Calculate total catch at equilibrium for all efforts
  Bfg.equil  <- colSums(B.ijte[,,t,]) # Calculate biomass for each functional group at equilibrium for all efforts
  cBfg.equil <- colSums(cB.ijte[,,t,]) # Calculate catch for each functional group at equilibrium for all efforts
  # Calculate mean weight. length and length at maturity
  BmeanW    <- vector(mode="numeric", length=length(effort))
  BmeanL    <- vector(mode="numeric", length=length(effort))
  BmeanLmat <- vector(mode="numeric", length=length(effort))
  CmeanW    <- vector(mode="numeric", length=length(effort))
  CmeanL    <- vector(mode="numeric", length=length(effort))
  for(e in 1:length(effort)){
    totalB       <- sum(B.ijte[,,t,e])
    totalN       <- sum(N.ijte[,,t,e])
    totalL       <- sum(rowSums(N.ijte[,,t,e]) * L.mid)
    BmeanW[e]    <- totalB/totalN
    BmeanL[e]    <- totalL/totalN
    BmeanLmat[e] <- sum(colSums(N.ijte[,,t,e]) * Lmat) / totalN
    totalcB      <- sum(cB.ijte[,,t,e])
    totalcN      <- sum(cN.ijte[,,t,e])
    totalcL      <- sum(rowSums(cN.ijte[,,t,e]) * L.mid)
    CmeanW[e]    <- totalcB/totalcN
    CmeanL[e]    <- totalcL/totalcN
  }
  # Calculate functional group specific mean weight, length, and length at maturity
  fg.BmeanW    <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.BmeanL    <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.BmeanLmat <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.CmeanW    <- array(NA, dim=c(1,nspecies,length(effort)))
  fg.CmeanL    <- array(NA, dim=c(1,nspecies,length(effort)))
  for(e in 1:length(effort)){
    fg.totalB <- colSums(B.ijte[,,t,e])
    fg.totalN <- colSums(N.ijte[,,t,e])
    fg.totalL <- colSums(N.ijte[,,t,e] * repmat(as.array(L.mid),1,nspecies))
    fg.BmeanW[,,e]    <- fg.totalB / fg.totalN
    fg.BmeanL[,,e]    <- fg.totalL / fg.totalN
    fg.BmeanLmat[,,e] <- (colSums(N.ijte[,,t,e]) * Lmat) / fg.totalN
    fg.totalcB <- colSums(cB.ijte[,,t,e])
    fg.totalcN <- colSums(cN.ijte[,,t,e])
    fg.totalcL <- colSums(cN.ijte[,,t,e] * repmat(as.array(L.mid),1,nspecies))
    fg.CmeanW[,,e] <- fg.totalcB/fg.totalcN
    fg.CmeanL[,,e] <- fg.totalcL/fg.totalcN
  }
  out <- list(B.equil,cB.equil,Bfg.equil,cBfg.equil,BmeanW,BmeanL,BmeanLmat,CmeanW,CmeanL,fg.BmeanW,fg.BmeanL,fg.BmeanLmat,fg.CmeanW,fg.CmeanL)
}

