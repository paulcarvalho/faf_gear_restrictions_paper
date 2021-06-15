# --------------------------------------------------------------------------------------------------------------------------
## File 'Main_Script.R'
## Gear-based fisheries population model
##
## Authors: Paul Carvalho
##          Austin Humphries
##
## Description: Coral reef fisheries model of nine functional groups to test biomass
##              and catch response to various gear-based management scenarios. The 
##              model losely represents the Wakatobi coral reef fishery in Indonesia,
##              where certain model parameters were derived.
# --------------------------------------------------------------------------------------------------------------------------

# Clean workspace
rm(list = ls())

# --------------------------------------------------- GENERAL NOTES --------------------------------------------------------
# 1. Use docstring([insert function name]) to view function documentation and information.

# --------------------------------------------------- LIBRARIES and DATA ---------------------------------------------------
library(readxl)
library(pracma)
library(ramify)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggstance)
library(ggpubr)
library(fishualize)
library(forcats)
library(fitdistrplus)
library(plotrix)
library(tidyverse)
library(splitstackshape)
library(data.table)
library(docstring)
library(roxygen2)
library(Rmisc)
library(foreach)
library(doParallel)

# Load data 
landings  <- read.csv("landings.csv") # fisheries-dependent data
uvc       <- read.csv("uvc.csv") # fisheries-independent data
LH.params <- read_xlsx("model_params.xlsx", sheet = 1) # life history parameters
foodweb   <- read_xlsx("model_params.xlsx", sheet = 2) # food web matrix

# --------------------------------------------------- MODEL PARAMETERS -----------------------------------------------------
t        <- 200 # number of timesteps
nspecies <- 9 # number of species/functional groups
nsc      <- 11 # number of size classes

# life-history parameters
k    <- LH.params$k # instantaneous growth parameters of the von Bertalanffy growth equation
Linf <- LH.params$Linf # asymptotic length parameters of the von Bertalanffy growth equation
W.a  <- LH.params$W_a # intercept parameters of length-weight relationship
W.b  <- LH.params$W_b # slope parameters of length-weight relationship
Lmin <- LH.params$Lmin # minimum sizes in centimeters
Lmax <- LH.params$Lmax # maximum sizes in centimeters
Lmat <- LH.params$Lmat # length at maturity

# lower and upper limits of size classes
maxsize <- max(Linf) # max Linf
minsize <- min(Lmin) # min fish size
L.lower <- seq(from = minsize, to = maxsize - (maxsize-minsize)/nsc, length.out = nsc)   # lower limit (included in the size class)
L.upper <- seq(from = minsize + ((maxsize-minsize)/nsc), to = maxsize, length.out = nsc) # upper limit (excluded)
L.mid   <- L.lower + (L.upper - L.lower)/2 # mid-point of each size class

# Beverton-Holt
h <- Lmat/Linf # dimensionless steepness parameter for Beverton-Holt stock-recruitment relationship

# Food web matrix
tau     <- as.matrix(foodweb[,c(2:10)]) # Foodweb matrix. Columns indicate predator species, rows indicate prey species, and 1 indicates prey consumed by predators
other.i <- 5e7

# Parameters for predator size preference function
mu    <- -3.5 # mean for lognormal function
sigma <- 1    # sd for lognormal function

# --------------------------------------------------- MODEL IMPLEMENTATION ---------------------------------------------------

source('model_functions.R')

# determine proportion of fish that grow to the next size class
phi.out <- calc_phi(L.lower, L.upper, Linf, k, nsc) 
phi     <- phi.out[[1]] # prop of fish that grow to next size class
phi.min <- phi.out[[2]] # length of a timestep (years)

# calculate ration
ration_out <- calc_ration(k, Linf, nsc, nspecies, L.lower, L.upper, L.mid, W.a, W.b, phi.min, scale.Ge = 0)
ration     <- ration_out[[1]]
weight     <- ration_out[[2]]
sc_Linf    <- ration_out[[3]]

# calculate size ratio based on preference parameters
M2_prefs <- calc_prefs(L.mid, nsc, nspecies, mu, sigma, weight, sc_Linf)

# calculate suitabilities for pred and prey
suit <- calc_suit(M2_prefs, tau, nsc, nspecies, sc_Linf)

# natural mortality rate (M1)
M1 <- nat_mortality(L.lower, L.upper, nspecies, nsc, phi.min, Linf, k, "mid") # natural mortality (excluding predation)

# bootstrap (resample) landings and uvc data to calculate base catchability*selectivity 
n.bs    <- 30
n.cores <- 5
cluster <- parallel::makeCluster(n.cores, type = "PSOCK"); registerDoParallel(cluster)
# ptm     <- proc.time()
q.tmp   <- foreach(i = 1:n.bs, .combine = 'rbind') %dopar% {
               library(dplyr)
               library(data.table)
               library(splitstackshape)
               library(plotrix)
               library(fitdistrplus)
               library(tidyr)
               q.tmp <- bootstrap_qs(landings, uvc, resample = TRUE)
               q.tmp <- array(as.numeric(unlist(q.tmp)), dim = c(nsc, nspecies, 3))
}
# proc.time() - ptm
parallel::stopCluster(cluster)
q.tmp2 <- array(NA, dim = c(nsc, nspecies, 3, n.bs))
for(i in 1:n.bs){
     tmp <- array(as.array(q.tmp)[i, ], dim = c(nsc, nspecies, 3))
     q.tmp2[, , , i] <- tmp 
}
test.hook  <- apply(q.tmp2[, , 1, ], c(1, 2), mean); test.hook  <- test.hook/sum(test.hook)
test.net   <- apply(q.tmp2[, , 2, ], c(1, 2), mean); test.net   <- test.net/sum(test.net)
test.spear <- apply(q.tmp2[, , 3, ], c(1, 2), mean); test.spear <- test.spear/sum(test.spear)
q <- array(c(test.hook, test.net, test.spear), dim = c(nsc, nspecies, 3)) 

# run fisheries model
source("run_model.R")

# --------------------------------------------------- CALIBRATE RECRUITMENT ---------------------------------------------------

# order: browser, detritivore, excavator/scraper, grazer, macro-invertivore, micro-invertivore, pisci-invertivore, piscivore, planktivore
bio.exp <- c(47, 87, 47, 47, 102, 102, 45, 45, 1067) # expected biomass/hectare in the absence of fishing; value calculated from Campbell et al. (2020) and within range reported in MacNeil et al. (2015)
r.i     <- c(169, 1525, 232, 284, 2211, 995, 35, 70, 6832) # initial estimate for recruits (values set equal to optimization)
Ntest   <- matrix(1, nsc, nspecies)
par.in  <- c(r.i, other.i) # optimize recruitment and food 'other' to achieve expected biomass of each functional group in the absence of fishing.

# ptm<-proc.time()
# r.fit <- optim(par=par.in, fn=cal_recruitment, N0=Ntest, t=t, nsc=nsc, nspecies=nspecies,
#             	 M1=M1, phi=phi, L.lower=L.lower, L.upper=L.upper, W.a=W.a, W.b=W.b, Bi.F0=bio.exp,
# 						   suit=suit, ration=ration, weight=weight, sc_Linf=sc_Linf, phi.min=phi.min,
#             	 control = list(maxit = 8000))
# proc.time()-ptm
# r <- r.fit$par[1:9] # c(169.92188, 1525.68866, 232.19953, 284.86446, 2211.87258, 995.01757, 35.14660, 70.13579, 6831.70676)
# other <- r.fit$par[10] # 55003764
r          <- c(169.92188, 1525.68866, 232.19953, 284.86446, 2211.87258, 995.01757, 35.14660, 70.13579, 6831.70676)
other      <- 55003764
alpha_beta <- r # asymptote of BH srr

# run model with calibrated recruitment and no fishing
N0         <- calc_N0(uvc)
N.ijt      <- array(NA, dim = c(nsc, nspecies, t))
N.ijt[,,1] <- N0
S.ijt      <- array(NA, dim = c(1, nspecies, t))
for(i in 2:t){
        N.tmp      <- N.ijt[,,i-1] # save temporary abundance from previous time step
        S.tmp      <- calc_ssb(nspecies, nsc, L.lower, L.upper, N.ijt, i, Lmat) # calculate spawning stock biomass
        N.tmp[1,]  <- N.tmp[1,] + r # implement recruitment
        M2         <- calc_M2(N=N.tmp, suit, ration, nspecies, nsc, other, weight, sc_Linf, phi.min) # calculate predation mortality
        N.tmp      <- N.tmp * exp(-(M1+M2)) # natural mortality
        N.tmp      <- calc_growth(N.tmp, phi, nsc, nspecies) # implement fish growth
        N.ijt[,,i] <- N.tmp # save abundance for current time step
        S.ijt[,,i] <- S.tmp # save spawning stock biomass for current time step
}

# Parameterize Beverton-Holt model 
SSB       <- calc_bio(S.ijt, t, L.lower, L.upper, W.a, W.b, 1, nspecies) # calculate SSB
Ro        <- r # estimated recruitment at virgin biomass
Bo        <- SSB[,,t]                     # estimated virgin spawning stock biomass
alpha     <- (Bo/Ro) * ((1-h)/4*h) # Beverton-Holt parameter (1/alpha is the maximum per capita production of recruits)
beta      <- (5*h - 1) / (4*h*Ro)   # Beverton-Holt parameter (R approaches 1/beta as biomass increases)
BH.df     <- data.frame(fg = rep(c("Browser", "Detritivore", "Excavator/Scraper", "Grazer", "Macro-invertivore", "Micro-invertivore", "Pisci-invertivore", "Piscivore", "Planktivore"), each = 1000), BHb = rep(seq(1,1000,1), 9), alpha = rep(alpha, each = 1000), beta = rep(beta, each = 1000))
BH.df$BHr <- BH.df$BHb / (BH.df$alpha + BH.df$beta * BH.df$BHb)
ggplot() + 
     geom_line(data=BH.df, aes(x = (BHb), y = (BHr), color = fg), size = 1, alpha = 0.8) + 
     labs(x = "Spawning stock biomass", y = "Recruits") + 
     scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")) + 
     scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + 
     theme_classic() + 
     theme(legend.title = element_blank(), legend.position = c(0.9,0.65))

# # get the abundance and biomass for each functional group in the absence of fishing
# Ni.F0 <- colSums(N.ijt[,,1]) # abundance
# N.it  <- as.data.frame(t(colSums(N.ijt)))
# N.it  <- N.it %>% mutate(time = 1:t) %>% dplyr::select(time, "Browser"=V1,"Detritivore"=V2,"Excavator/scraper"=V3,"Grazer"=V4,"Macro-invertivore"=V5, "Micro-invertivore"=V6,"Pisci-invertivore"=V7,"Piscivore"=V8,"Planktivore"=V9) %>% gather(., "fg", "N", -time)
# # calculate biomass
# B.ijt <- calc_bio(N.ijt, t, L.lower, L.upper, W.a, W.b, nsc, nspecies)
# Bi.F0 <- colSums(B.ijt[,,t]) # biomass
# B.it  <- as.data.frame(t(colSums(B.ijt)))
# B.it  <- B.it %>% mutate(time = 1:t) %>% dplyr::select(time, "Browser"=V1,"Detritivore"=V2,"Excavator/scraper"=V3,"Grazer"=V4,"Macro-invertivore"=V5, "Micro-invertivore"=V6,"Pisci-invertivore"=V7,"Piscivore"=V8,"Planktivore"=V9) %>% gather(., "fg", "B", -time)
# # plot 
# plot.equil1 <- ggplot() + geom_line(data = N.it, aes(x = time, y = log(N), color = fg), lwd = 1) + labs(y = "log(Abundance)") + scale_y_continuous(limits = c(min(log(N.it$N)),9), expand = c(0,0)) + theme_classic() + guides(color=guide_legend(ncol=3)) + theme(legend.title = element_blank(), legend.position = c(0.5,0.9))
# plot.equil2 <- ggplot() + geom_line(data = B.it, aes(x = time, y = log(B), color = fg), lwd = 1) + labs(x = "Timesteps", y = "log(Biomass)") + theme_classic() + theme(legend.position = "none")
# ggarrange(plot.equil1, plot.equil2, ncol = 1, nrow = 2, labels = c("A", "B"))

# --------------------------------------------------- MODEL SETTINGS ---------------------------------------------------

max.effort <- 29 # maximum effort for each gear type. Units are arbitrary
effort     <- seq(0, 1, length.out = 20) * max.effort # fishing effort to test
nbs        <- 10 # number of times to run bootstrap

# --------------------------------------------------- MODEL SCENARIO 1: ALL GEARS ---------------------------------------------------
gear.mgmt.1 <- c(1,1,1) # all gears used
scen.1      <- run_model(effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min, N0)
scen.1.out  <- calc_summary_indices(N.ijte = scen.1[[1]], B.ijte = scen.1[[2]], cN.ijte = scen.1[[3]], cB.ijte = scen.1[[4]], t, L.mid, Lmat, nspecies)
effort.1.bs <- c(7.25, 14.5, 21.75)

scen.1.bs   <- run_bsmodel(nbs, landings, uvc, effort.1.bs, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min, N0)
plot1 <- ggplot() + 
     geom_line(aes(x = effort/max.effort, y = scen.1.out[[1]])) +
     geom_point(aes(x=effort.1.bs/max.effort, y=scen.1.bs[[1]]$B.mu)) +
     geom_errorbar(aes(x=effort.1.bs/max.effort,
                       ymin=scen.1.bs[[1]]$B.lo,
                       ymax=scen.1.bs[[1]]$B.up), width = 0.05)

# --------------------------------------------------- MODEL SCENARIO 2: NO LINE FISHING ---------------------------------------------------
gear.mgmt.2 <- c(0,1,1) # No hook-and-line fishing
scen.2      <- run_model(effort, gear.mgmt.2, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min)
scen.2.out  <- calc_summary_indices(N.ijte = scen.2[[1]], B.ijte = scen.2[[2]], cN.ijte = scen.2[[3]], cB.ijte = scen.2[[4]],  t, L.mid, Lmat, nspecies)
effort.2.bs <- c(7.25, 14.5, 21.75)

scen.2.bs   <- run_bsmodel(nbs, landings, uvc, effort.2.bs, gear.mgmt.2, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min)
plot1 +
     geom_line(aes(x = effort/max.effort, y = scen.2.out[[1]]), color='red') +
     geom_point(aes(x = effort.2.bs/max.effort, y = scen.2.bs[[1]]$B.mu), color='red') +
     geom_errorbar(aes(x=effort.2.bs/max.effort,
                       ymin=scen.2.bs[[1]]$B.lo,
                       ymax=scen.2.bs[[1]]$B.up), width = 0.05, color='red')


# --------------------------------------------------- MODEL SCENARIO 3: NO NET FISHING ---------------------------------------------------
gear.mgmt.3 <- c(1,0,1) # No net fishing
scen.3      <- run_model(effort, gear.mgmt.3, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min)
scen.3.out  <- calc_summary_indices(N.ijte = scen.3[[1]], B.ijte = scen.3[[2]], cN.ijte = scen.3[[3]], cB.ijte = scen.3[[4]], t, L.mid, Lmat, nspecies)
effort.3.bs <- c(7.25, 14.5, 21.75)

# --------------------------------------------------- MODEL SCENARIO 4: NO SPEAR FISHING ---------------------------------------------------
gear.mgmt.4 <- c(1,1,0) # No spear fishing
scen.4      <- run_model(effort, gear.mgmt.4, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min)
scen.4.out  <- calc_summary_indices(N.ijte = scen.4[[1]], B.ijte = scen.4[[2]], cN.ijte = scen.4[[3]], cB.ijte = scen.4[[4]], t, L.mid, Lmat, nspecies)
effort.4.bs <- c(7.25, 14.5, 21.75)

# --------------------------------------------------- MODEL SCENARIO 5: LINE FISHING ---------------------------------------------------
gear.mgmt.5 <- c(1,0,0) # Only hook-and-line fishing
scen.5      <- run_model(effort, gear.mgmt.5, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min)
scen.5.out  <- calc_summary_indices(N.ijte = scen.5[[1]], B.ijte = scen.5[[2]], cN.ijte = scen.5[[3]], cB.ijte = scen.5[[4]], t, L.mid, Lmat, nspecies)
effort.5.bs <- c(7.25, 14.5, 21.75)

# --------------------------------------------------- MODEL SCENARIO 6: NET FISHING ---------------------------------------------------
gear.mgmt.6 <- c(0,1,0) # Only net fishing
scen.6      <- run_model(effort, gear.mgmt.6, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min)
scen.6.out  <- calc_summary_indices(N.ijte = scen.6[[1]], B.ijte = scen.6[[2]], cN.ijte = scen.6[[3]], cB.ijte = scen.6[[4]],  t, L.mid, Lmat, nspecies)
effort.6.bs <- c(7.25, 14.5, 21.75)

# --------------------------------------------------- MODEL SCENARIO 7: SPEAR FISHING ---------------------------------------------------
gear.mgmt.7 <- c(0,0,1) # Only spear fishing
scen.7      <- run_model(effort, gear.mgmt.7, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta, suit, ration, other, weight, sc_Linf, phi.min)
scen.7.out  <- calc_summary_indices(N.ijte = scen.7[[1]], B.ijte = scen.7[[2]], cN.ijte = scen.7[[3]], cB.ijte = scen.7[[4]], t, L.mid, Lmat, nspecies)
effort.7.bs <- c(7.25, 14.5, 21.75)

# --------------------------------------------------- MODEL SENSITIVITY ANALYSES ---------------------------------------------------
# First set of sensitivity analyses run simulations with key parameters (i.e., alpha, beta, mu, sigma, Ge, and tau) increased or
# decreased by 10%. Then indicators for each sensitivity analysis are calculated with cal_sensitivity_indices().

# Second, sensitivity to gear specifications is tested. We made the simplifying assumption that gear size (i.e., hook gape size and
# net mesh size) is proportional to the lognormal mean of size selectivity. We then calculated the lognormal means that indicate
# 1 and 2 cm increases in hook gape size and net mesh size. Simulations for modified gear specifications were only run with all 
# fishing gears used simultaneously.

# --------------------------------------------------- MODEL SENSITIVITY - NO FISHING ---------------------------------------------------

base.effort <- effort

base.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                        suit, ration, other, weight, sc_Linf, phi.min)


base.sa <- calc_sensitivity_indices(N.ijte = base.model[[1]], B.ijte = base.model[[2]], 
							 cN.ijte = base.model[[3]], cB.ijte = base.model[[4]], 
							 L.mid, Lmat, tau, base.effort)


# --------------------------------------------------- MODEL SENSITIVITY - alpha ---------------------------------------------------

alpha.hi <- alpha + alpha*0.1
alpha.lo <- alpha - alpha*0.1

alphaHi.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha.hi, beta,
                           suit, ration, other, weight, sc_Linf, phi.min)

alphaHi.sa <- calc_sensitivity_indices(N.ijte = alphaHi.model[[1]], B.ijte = alphaHi.model[[2]], 
							    cN.ijte = alphaHi.model[[3]], cB.ijte = alphaHi.model[[4]],
							    L.mid, Lmat, tau, base.effort)

alphaLo.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha.lo, beta,
                           suit, ration, other, weight, sc_Linf, phi.min)

alphaLo.sa <- calc_sensitivity_indices(N.ijte = alphaLo.model[[1]], B.ijte = alphaLo.model[[2]], 
							    cN.ijte = alphaLo.model[[3]], cB.ijte = alphaLo.model[[4]],
							    L.mid, Lmat, tau, base.effort)


# --------------------------------------------------- MODEL SENSITIVITY - beta ---------------------------------------------------

beta.hi <- beta + beta*0.1
beta.lo <- beta - beta*0.1

betaHi.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta.hi,
                          suit, ration, other, weight, sc_Linf, phi.min)

betaHi.sa <- calc_sensitivity_indices(N.ijte = betaHi.model[[1]], B.ijte = betaHi.model[[2]], 
							   cN.ijte = betaHi.model[[3]], cB.ijte = betaHi.model[[4]],
							   L.mid, Lmat, tau, base.effort)

betaLo.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta.lo,
                          suit, ration, other, weight, sc_Linf, phi.min)

betaLo.sa <- calc_sensitivity_indices(N.ijte = betaLo.model[[1]], B.ijte = betaLo.model[[2]], 
							   cN.ijte = betaLo.model[[3]], cB.ijte = betaLo.model[[4]],
							   L.mid, Lmat, tau, base.effort)

# --------------------------------------------------- MODEL SENSITIVITY - mu ---------------------------------------------------

mu.hi <- mu + mu*0.1
mu.lo <- mu - mu*0.1

M2_prefs_muHi <- calc_prefs(L.mid, nsc, nspecies, mu.hi, sigma, weight, sc_Linf)
M2_prefs_muLo <- calc_prefs(L.mid, nsc, nspecies, mu.lo, sigma, weight, sc_Linf)

suit_muHi <- calc_suit(M2_prefs_muHi, tau, nsc, nspecies, sc_Linf)
suit_muLo <- calc_suit(M2_prefs_muLo, tau, nsc, nspecies, sc_Linf)

muHi.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                        suit_muHi, ration, other, weight, sc_Linf, phi.min)

muHi.sa <- calc_sensitivity_indices(N.ijte = muHi.model[[1]], B.ijte = muHi.model[[2]], 
							 cN.ijte = muHi.model[[3]], cB.ijte = muHi.model[[4]],
							 L.mid, Lmat, tau, base.effort)

muLo.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                        suit_muLo, ration, other, weight, sc_Linf, phi.min)

muLo.sa <- calc_sensitivity_indices(N.ijte = muLo.model[[1]], B.ijte = muLo.model[[2]], 
							 cN.ijte = muLo.model[[3]], cB.ijte = muLo.model[[4]],
							 L.mid, Lmat, tau, base.effort)


# --------------------------------------------------- MODEL SENSITIVITY - sigma ---------------------------------------------------

sigma.hi <- sigma + sigma*0.1
sigma.lo <- sigma - sigma*0.1

M2_prefs_sigmaHi <- calc_prefs(L.mid, nsc, nspecies, mu, sigma.hi, weight, sc_Linf)
M2_prefs_sigmaLo <- calc_prefs(L.mid, nsc, nspecies, mu, sigma.lo, weight, sc_Linf)

suit_sigmaHi <- calc_suit(M2_prefs_sigmaHi, tau, nsc, nspecies, sc_Linf)
suit_sigmaLo <- calc_suit(M2_prefs_sigmaLo, tau, nsc, nspecies, sc_Linf)

sigmaHi.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                           suit_sigmaHi, ration, other, weight, sc_Linf, phi.min)

sigmaHi.sa <- calc_sensitivity_indices(N.ijte = sigmaHi.model[[1]], B.ijte = sigmaHi.model[[2]], 
							    cN.ijte = sigmaHi.model[[3]], cB.ijte = sigmaHi.model[[4]],
							    L.mid, Lmat, tau, base.effort)

sigmaLo.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                           suit_sigmaLo, ration, other, weight, sc_Linf, phi.min)

sigmaLo.sa <- calc_sensitivity_indices(N.ijte = sigmaLo.model[[1]], B.ijte = sigmaLo.model[[2]], 
							    cN.ijte = sigmaLo.model[[3]], cB.ijte = sigmaLo.model[[4]],
							    L.mid, Lmat, tau, base.effort)


# --------------------------------------------------- MODEL SENSITIVITY - Ge ---------------------------------------------------

ration_outGeHi <- calc_ration(k, Linf, nsc, nspecies, L.lower, L.upper, L.mid, W.a, W.b, phi.min, scale.Ge = 0.1)
ration.GeHi <- ration_outGeHi[[1]]

ration_outGeLo <- calc_ration(k, Linf, nsc, nspecies, L.lower, L.upper, L.mid, W.a, W.b, phi.min, scale.Ge = -0.1)
ration.GeLo <- ration_outGeLo[[1]]

GeHi.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                        suit, ration.GeHi, other, weight, sc_Linf, phi.min)

GeHi.sa <- calc_sensitivity_indices(N.ijte = GeHi.model[[1]], B.ijte = GeHi.model[[2]], 
							 cN.ijte = GeHi.model[[3]], cB.ijte = GeHi.model[[4]],
							 L.mid, Lmat, tau, base.effort)

GeLo.model <- run_model(base.effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q, alpha, beta,
                        suit, ration.GeLo, other, weight, sc_Linf, phi.min)

GeLo.sa <- calc_sensitivity_indices(N.ijte = GeLo.model[[1]], B.ijte = GeLo.model[[2]], 
							 cN.ijte = GeLo.model[[3]], cB.ijte = GeLo.model[[4]],
							 L.mid, Lmat, tau, base.effort)


# --------------------------------------------------- MODEL SENSITIVITY - gear specification ---------------------------------------------------

# Increase mean of the lognormal function for size selectivity of hook-and-line and net fishing.
# Simulating an increase in hook and mesh size.

# Hook sizes
hooksize.i <- 1.5 # cm of hook gape size
hooksize.1 <- hooksize.i + 1 # 1 cm increase in hook gape size
hooksize.2 <- hooksize.i + 2 # 2 cm increase in hook gape size

# Mesh sizes
netsize.i <- 6.5 # cm of mesh diagonal
netsize.1 <- netsize.i + 1 # 1 cm increase in mesh size
netsize.2 <- netsize.i + 2 # 2 cm increase in mesh size

# Selectivity for a 1 cm increase in hook and net mesh size
mu.hook1 <- log((exp(size_sel$mu[1]) / hooksize.i) * hooksize.1)
mu.net1 <- log((exp(size_sel$mu[2]) / netsize.i) * netsize.1)

# Selectivity for a 2 cm increase in hook and net mesh size
mu.hook2 <- log((exp(size_sel$mu[1]) / hooksize.i) * hooksize.2)
mu.net2 <- log((exp(size_sel$mu[2]) / netsize.i) * netsize.2)

# Calculate catchability and selectivity for 1 cm increase
q1.sa <- array(NA, dim = c(nsc, nspecies, 3)) # create 3D array with size classes (rows), species (functional groups; cols), and gear types.
q1.sa[,,1] <- calc_qs(L.mid, mu.hook1, size_sel$sd[1], fg_cat$q_line) # INCREASED hook size
q1.sa[,,2] <- calc_qs(L.mid, mu.net1, size_sel$sd[2], fg_cat$q_net) # INCREASED net size
q1.sa[,,3] <- calc_qs(L.mid, size_sel$mu[3], size_sel$sd[3], fg_cat$q_spear) # spear

# Calculate catchability and selectivity for 2 cm increase
q2.sa <- array(NA, dim = c(nsc, nspecies, 3)) # create 3D array with size classes (rows), species (functional groups; cols), and gear types.
q2.sa[,,1] <- calc_qs(L.mid, mu.hook2, size_sel$sd[1], fg_cat$q_line) # INCREASED hook size
q2.sa[,,2] <- calc_qs(L.mid, mu.net2, size_sel$sd[2], fg_cat$q_net) # INCREASED net size
q2.sa[,,3] <- calc_qs(L.mid, size_sel$mu[3], size_sel$sd[3], fg_cat$q_spear) # spear

# MODEL WITH 1 CM INCREASE
# gear.mgmt.1 <- c(1,1,1) # all gears used
q1sa.model <- run_model(effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q1.sa, alpha, beta,
					    					suit, ration, other, weight, sc_Linf, phi.min)

q1sa.out <- calc_summary_indices(N.ijte = q1sa.model[[1]], B.ijte = q1sa.model[[2]], cN.ijte = q1sa.model[[3]], cB.ijte = q1sa.model[[4]], 
                                 t, L.mid, Lmat, nspecies)

# MODEL WITH 2 CM INCREASE
# gear.mgmt.1 <- c(1,1,1) # all gears used
q2sa.model <- run_model(effort, gear.mgmt.1, nsc, nspecies, t, Lmat, M1, phi, L.lower, L.upper, W.a, W.b, q2.sa, alpha, beta,
					    					suit, ration, other, weight, sc_Linf, phi.min)

q2sa.out <- calc_summary_indices(N.ijte = q2sa.model[[1]], B.ijte = q2sa.model[[2]], cN.ijte = q2sa.model[[3]], cB.ijte = q2sa.model[[4]], 
                                 t, L.mid, Lmat, nspecies)

# --------------------------------------------------- SAVE DATA ---------------------------------------------------

save.image(file = "model_data.RData")


