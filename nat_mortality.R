### nat_mortality
###
### Author: Paul Carvalho
###
### Description: Calculates natural mortality rates for each species and size class
###              based on equation from Gislason et al. (2010)

# ln(Mi,j) = 0.55 - 1.61 ln(Li,j) + 1.44 ln(Linf) + ln(ki)

nat_mortality <- function(L.lower, L.upper, nspecies, nsc, phi.min, Linf, k, which.L){
        # create matrices for VBGF params
        Linf.mat <- repmat(t(as.matrix(Linf)), nsc, 1)
        k.mat <- repmat(t(as.matrix(k)), nsc, 1)
        
        if(which.L == "mid"){
                # calculate midpoint of each size class
                Lmid <- (L.lower + L.upper) / 2
                Lmid.mat <- repmat(as.matrix(Lmid), 1, nspecies)
    
                # calculate natural mortality and scale to time step of fastest growing species-size class combination
                lnM <- 0.55 - (1.61*log(Lmid.mat)) + (1.44*log(Linf.mat)) + (log(k.mat))
                M.year <- exp(lnM)
                M <- M.year*phi.min # scale to timestep
                
                # Plot natural mortality curves
                cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
                M.df <- as.data.frame(M)
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
        } else if(which.L == "upper") {
                # create matrix
                Lup.mat <- repmat(as.matrix(L.upper), 1, nspecies)
    
                # calculate natural mortality and scale to time step of fastest growing species-size class combination
                lnM <- 0.55 - (1.61*log(Lup.mat)) + (1.44*log(Linf.mat)) + (log(k.mat))
                M.year <- exp(lnM)
                M <- M.year*phi.min # scale to timestep
        } else if(which.L == "lower") {
                # create matrix
                Llo.mat <- repmat(as.matrix(L.lower), 1, nspecies)
     
                # calculate natural mortality and scale to time step of fastest growing species-size class combination
                lnM <- 0.55 - (1.61*log(Llo.mat)) + (1.44*log(Linf.mat)) + (log(k.mat))
                M.year <- exp(lnM)
                M <- M.year*phi.min # scale to timestep
        } else if(which.L == "mean") {
                L.add <- 0:((L.upper[1] - L.lower[1]) - 1)
                M.tmp <- array(NA, dim = c(nsc, nspecies, length(L.add)))
                
                for(i in 1:length(L.add)){
                    L <- L.lower + L.add[i]
                    L <- repmat(as.matrix(L), 1, nspecies)
                    lnM <- 0.55 - (1.61*log(L)) + (1.44*log(Linf.mat)) + (log(k.mat))
                    M.year <- exp(lnM)
                    M.tmp1 <- M.year*phi.min # scale to timestep
                    M.tmp[,,i] <- M.tmp1
                }
                M <- apply(M.tmp, c(1,2), mean)
        }
        
        return(M)
}
