### size_dist
###
### Author: Paul Carvalho
###
### Description: This script is called by "main_model.R" and generates the size distribution
###              using the inverse method of the power-law density function in 
###              Edwards et al. (2017).


# Libraries ---------------------------------------------------------------
library(DescTools)

# Functions ---------------------------------------------------------------
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

# Length spectra ----------------------------------------------------------
# Length spectra values. Values fall within estimates derived in Robinson et al. (2017).
# Herbivores and detritivores were set shallower than other functional groups due to
# metabolic constraints on abundance-size relationship (Trebilco et al. 2013,
# Robinson and Baum 2019).
# b = -2.3 # herbivores and detritivores
# b = -2.7 # piscivores and invertivores and planktivores

fg.abundance <- as.data.frame(mean.abundance.env$fg.abundance) # from mean_abundance.R

fg <- c("Browser","Detritivore","Excavator/scraper","Grazer","Macro-invertivore",
        "Micro-invertivore","Pisci-invertivore","Piscivore","Planktivore")

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

