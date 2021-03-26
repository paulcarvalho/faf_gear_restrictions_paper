### bootstrap_qs
###
### Author: Paul Carvalho
###
### Description: Re-sample (via bootstrapping) landings and UVC data to quantify uncertainty in selectivity values

bootstrap_qs <- function(landings, uvc, resample = TRUE){
  # Filter out size that are not included in the model and not exploited in the fishery
  landings <- landings %>% filter(size_cm >= 10 & size_cm < 65)
  uvc      <- uvc %>% filter(size_cm >= 10 & size_cm < 65)
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
      tmp.bs$site.tran <- site.sample$site_tran[i] # save the site name and transect
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
  n_i         <- uvc %>% dplyr::group_by(site.tran, fg) %>% mutate(abundance = 1) %>% summarise(abundance = sum(abundance)) %>% dplyr::group_by(site.tran) %>% mutate(ratio = abundance/sum(abundance)) %>% dplyr::group_by(fg) %>% summarise(n = mean(ratio), n_stderr = std.error(ratio))
  names(n_i)  <- c("fg", "n", "n_stderr")
  n_ij        <- uvc %>% dplyr::group_by(site.tran, fg, size_5cm_bin) %>% mutate(abundance = 1) %>% summarise(abundance = sum(abundance)) %>% dplyr::group_by(site.tran) %>% mutate(ratio = abundance/sum(abundance)) %>% dplyr::group_by(fg, size_5cm_bin) %>% summarise(n = mean(ratio), n_stderr = std.error(ratio))
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



