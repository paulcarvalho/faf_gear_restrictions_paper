# Functions
# Author: Paul Carvalho


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
      mutate(size = Lmid) %>%
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

bootstrap_qs <- function(landings, uvc, resample = TRUE){
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
  landings <- landings %>% filter(size_cm >= 10 & size_cm < 65)
  uvc      <- uvc %>% filter(size_cm >= 10 & size_cm < 65) %>% mutate(site_tran = paste(site_name, transect))
  uni.fg   <- unique(landings$fg) # Get unique functional groups
  
  if(resample == TRUE){
    # Resample landings data via bootstrapping
    trip.sample <- landings %>% group_by(trip_id, gear_cat1) %>% summarize(size = sum(abundance)) # Save trip IDs and gear and number of fish caught
    landings.bs <- NULL # empty dataframe to fill during resampling of landings data
    for(i in 1:length(trip.sample$trip_id)){ # Iterate through trips and resample from fish caught by each gear
      gear.i           <- trip.sample$gear_cat1[i] # save gear for trip.id[i]
      sample.i         <- landings %>% filter(gear_cat1 == gear.i) %>% dplyr::select(fg, size_cm, bin_5cm) %>% filter(!(size_cm < 10)) # save all fg and sizes caught by gear.i
      index.bs         <- sample(x = 1:length(sample.i$fg), size = trip.sample$size[i], replace = TRUE) # sample indexes from all landings for gear.i
      tmp.bs           <- sample.i[index.bs, ] # save the functional groups and size classes that were sampled
      tmp.bs$trip_id   <- trip.sample$trip_id[i] # save trip
      tmp.bs$gear_cat1 <- trip.sample$gear_cat1[i] # save gear
      landings.bs      <- rbind(landings.bs, tmp.bs) # add to full bootstrapped sample
    }
    landings <- landings.bs
    
    # Resample UVC data via bootstrapping
    site.sample <- uvc %>% mutate(site_tran = paste(site_name, transect)) %>% dplyr::select(site_name, site_tran, fg, size_5cm_bin, abundance) %>% group_by(site_name, site_tran) %>% summarize(size = sum(abundance)) # Save site and transects and number of fish observed
    uvc.bs      <- NULL # empty dataframe to fill during resampling of landings data
    for(i in 1:length(site.sample$site_tran)){ # Iterate through site_trans and resample from fish observed at each site
      site.i           <- site.sample$site_name[i] # save site name for site_tran[i]
      sample.i         <- uvc %>% filter(site_name == site.i) %>% dplyr::select(fg, size_5cm_bin, abundance) # save all fish observed at site.i
      sample.i         <- setDT(expandRows(sample.i, 'abundance')) # expand table based on abundance
      index.bs         <- sample(x = 1:length(sample.i$fg), size = site.sample$size[i], replace = TRUE) # sample indexes from all fish observed at site.i
      tmp.bs           <- sample.i[index.bs, ] # save the functional groups and size classes that were sampled
      tmp.bs$site      <- site.i # save site name
      tmp.bs$site_tran <- site.sample$site_tran[i] # save the site name and transect
      uvc.bs           <- rbind(uvc.bs, tmp.bs) # add to full bootstrapped sample
    }
    uvc <- uvc.bs
  }
  
  # Calculate the ratio of catch per trip for each functional group (i), and functional group and size class (ij) - bootstrapped landings data
  c_i         <- landings %>% dplyr::group_by(trip_id, gear_cat1, fg) %>% mutate(abundance = 1) %>% summarise(abundance = sum(abundance)) %>% dplyr::group_by(trip_id, gear_cat1) %>% mutate(ratio = abundance/sum(abundance)) %>% dplyr::group_by(fg, gear_cat1) %>% summarise(c = mean(ratio), c_stderr = std.error(ratio))
  names(c_i)  <- c("fg", "gear", "c", "c_stderr")
  c_ij        <- landings %>% dplyr::group_by(trip_id, gear_cat1, fg, bin_5cm) %>% mutate(abundance = 1) %>% summarise(abundance = sum(abundance)) %>% dplyr::group_by(trip_id, gear_cat1) %>% mutate(ratio = abundance/sum(abundance)) %>% dplyr::group_by(fg, bin_5cm, gear_cat1) %>% summarise(c = mean(ratio), c_stderr = std.error(ratio))
  names(c_ij) <- c("fg", "sc", "gear", "c", "c_stderr")
  
  # Calculate the ratio of numbers per transect per functional group (i) and functional group and size class (ij) - bootstrapped UVC data  
  n_i         <- uvc %>% dplyr::group_by(site_tran, fg) %>% mutate(abundance = 1) %>% summarise(abundance = sum(abundance)) %>% dplyr::group_by(site_tran) %>% mutate(ratio = abundance/sum(abundance)) %>% dplyr::group_by(fg) %>% summarise(n = mean(ratio), n_stderr = std.error(ratio))
  names(n_i)  <- c("fg", "n", "n_stderr")
  n_ij        <- uvc %>% dplyr::group_by(site_tran, fg, size_5cm_bin) %>% mutate(abundance = 1) %>% summarise(abundance = sum(abundance)) %>% dplyr::group_by(site_tran) %>% mutate(ratio = abundance/sum(abundance)) %>% dplyr::group_by(fg, size_5cm_bin) %>% summarise(n = mean(ratio), n_stderr = std.error(ratio))
  names(n_ij) <- c("fg", "sc", "n", "n_stderr")  
  
  cn_ij <- c_ij %>% full_join(., n_ij, by = c("fg", "sc")) # Merge landings and UVC dataframes
  q_i   <- c_i %>% full_join(., n_i, by = c("fg")) %>% mutate(q = c/n) # Calculate catchability; fisheries-dependent/-independent
  L.lower <- seq(from = 0, to = 65 - (65-5)/12, length.out = 13)   # low limit included in the size class
  L.upper <- seq(from = 0 + ((65-5)/12), to = 65, length.out = 13) # upper limit excluded
  Lmid <- ((L.lower + L.upper) / 2)[2:12] # Only include exploited size classes or modeled size classes
  
  hook.dist  <- landings %>% filter(gear_cat1 == "line");  hf.lnorm <- fitdist(hook.dist$size_cm, "lnorm") # fit lognormal distribution to hook-and-line landings
  net.dist   <- landings %>% filter(gear_cat1 == "net");   nf.lnorm <- fitdist(net.dist$size_cm, "lnorm") # fit lognormal distribution to net landings data
  spear.dist <- landings %>% filter(gear_cat1 == "spear"); sf.lnorm <- fitdist(spear.dist$size_cm, "lnorm") # fit lognormal distribution to spear landings data
  
  # Optimization function
  q_optim <- function(params, x, mu, sd, cn_ij, gear_type, fg, sc){
    q <- params[1:9] # functional group catchabilities, calibrated through optimization
    # model
    sj    <- dlnorm(x, meanlog = mu, sdlog = sd)
    sj    <- as.matrix(sj)
    model <- sj %*% q
    # expected
    cn_ij_x    <- cn_ij %>% mutate(fg = as.factor(fg)) %>% filter(gear == gear_type)
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
  hook_par <- t(as.matrix(par$q)) # save catchabilities for optimization
  mu.hook  <- hf.lnorm$estimate[1] # log(mean) for size selectivity
  sd.hook  <- hf.lnorm$estimate[2] # log(standard deviation) for size selectivity
  x_hook   <- Lmid # Set size classes to model
  hook_fit <- optim(par = hook_par, fn = q_optim, lower = c(0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), upper = c(1, 500, 500, 500, 500, 500, 500, 500, 500), method = "L-BFGS-B", control = list(maxit = 5000), x = x_hook, mu = mu.hook, sd = sd.hook, cn_ij = cn_ij, gear_type = "line", fg = unique(landings$fg), sc = unique(landings$bin_5cm))
  F_hook             <- (as.matrix(dlnorm(x_hook,mu.hook,sd.hook))) %*% (as.matrix(hook_fit$par)) # Calculate q*s with optimized parameters
  F_hook[F_hook < 0] <- 0 # Set negative values, if any, to zero
  F_hook             <- F_hook / sum(F_hook) # Calculate proportional q*s values
  # hook_heat          <- heatmap(t(F_hook), Colv=NA, Rowv=NA, scale="none", labCol=levels(cn_ij$sc), labRow=levels(as.factor(cn_ij$fg)), margins=c(9.5,9.5))
  
  # Net catchability (q) optimization
  par     <- subset(q_i, gear == "net")[,c(1,7)] # initial catchabilities for net fishing
  missing <- setdiff(uni.fg, par$fg) # check if functional groups are missing
  if (length(missing) > 0) {par <- rbind(data.frame(fg = missing, q = 0), (as.data.frame(par)))} # fill in missing functional group(s) with q=0
  net_par <- t(as.matrix(par$q)) # save catchabilities for optimization
  mu.net  <- nf.lnorm$estimate[1] # log(mean) for size selectivity
  sd.net  <- nf.lnorm$estimate[2] # log(standard deviation) for size selectivity
  x_net   <- Lmid # Set size classes to model
  net_fit <- optim(par = net_par, fn = q_optim, lower = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), upper = c(500, 500, 500, 500, 500, 500, 500, 500, 500), method = "L-BFGS-B", control = list(maxit = 10000), x = x_net, mu = mu.net, sd = sd.net, cn_ij = cn_ij, gear_type = "net", fg = unique(landings$fg), sc = unique(landings$bin_5cm))
  F_net            <- (as.matrix(dlnorm(x_net, mu.net, sd.net))) %*% (as.matrix(net_fit$par)) # Calculate q*s with optimized parameters
  F_net[F_net < 0] <- 0 # Set negative values, in any, to zero
  F_net            <- F_net / sum(F_net) # Caclulate proportional q*s values
  # net_heat <- heatmap(t(F_net), Colv=NA, Rowv=NA, scale="none", labCol=levels(cn_ij$sc), labRow=levels(as.factor(cn_ij$fg)), margins=c(9.5,9.5))
  
  # Spear catchability (q) optimization
  par       <- subset(q_i, gear == "spear")[,c(1,7)] # initial catchabilities for spear fishing
  missing   <- setdiff(uni.fg, par$fg) # check if functional groups are missing
  if (length(missing) > 0) {par <- rbind(data.frame(fg = missing, q = 0), (as.data.frame(par)))} # fill in missing functional group(s) with q=0
  spear_par <- t(as.matrix(par$q)) # save catchabilities for optimization
  mu.spear  <- sf.lnorm$estimate[1] # log(mean) for size selectivity
  sd.spear  <- sf.lnorm$estimate[2] # log(standard deviation) for size selectivity
  x_spear   <- Lmid # Set size classes to model  
  spear_fit <- optim(par = spear_par, fn = q_optim, lower = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), upper = c(500, 500, 500, 500, 500, 500, 500, 500, 500), method = "L-BFGS-B", control = list(maxit = 10000), x = x_spear, mu = mu.spear, sd = sd.spear, cn_ij = cn_ij, gear_type = "spear", fg = unique(landings$fg), sc = unique(landings$bin_5cm))
  F_spear              <- (as.matrix(dlnorm(x_spear, mu.spear, sd.spear))) %*% (as.matrix(spear_fit$par)) # Calculate q*s with optimized parameters
  F_spear[F_spear < 0] <- 0 # Set negative values, if any, to zero
  F_spear              <- F_spear / sum(F_spear) # Calculate proportional q*s values
  # spear_heat <- heatmap(t(F_spear), Colv=NA, Rowv=NA, scale="none", labCol=levels(cn_ij$sc), labRow=levels(as.factor(cn_ij$fg)), margins=c(9.5,9.5))
  
  return(list(F_hook, F_net, F_spear))  
}
