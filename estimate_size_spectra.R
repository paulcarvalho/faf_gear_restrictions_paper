### estimate_size_spectra
###
### Author: Paul Carvalho
###
### Description: Calculate the size spectra for each functional group using
###              methods outlined in Edwards et al. (2017).

# Libraries ---------------------------------------------------------------
# devtools::install_github("andrew-edwards/sizeSpectra")
library(sizeSpectra)
library(data.table)
library(splitstackshape)

# Expand dataframe such that abundance is 1 for each row
df <- setDT(expandRows(uvc.fish,"abundance"))[, abundance := sprintf("1")][]
uvc.fish.1 <- as.data.frame(df)

# Convert abundance to biomass
uvc.fish.2 <- uvc.fish.1 %>%
     mutate(abundance = as.numeric(abundance)) %>%
     mutate(biomass_g = ((a * (size_cm ^ b)) * as.numeric(abundance))) %>%
     mutate(biomass_kg = round((biomass_g/1000), digits = 2)) %>%
     filter(size_cm > 5 & size_cm < 67) %>% # underwater visual census is best for fishes in this size range (Kulbicki 1998)
     filter(family != "caesionidae") %>% # caesionids found in large school and observers have difficultly estimating abundance (Samoilys and Carlos 2000)
     filter(biomass_kg != 0) %>%
     filter(biomass_kg > 0.03) %>% # fish smaller than this size are likely inadequately sampled
     filter(region == "wakatobi")

# subset data by functional group, representative species
pisc <- uvc.fish.2 %>% filter(func_group == "Piscivore")
inve <- uvc.fish.2 %>% filter(func_group == "Invertivore")
herb <- uvc.fish.2 %>% filter(func_group == "Herbivore")
detr <- uvc.fish.2 %>% filter(func_group == "Detritivore")
plan <- uvc.fish.2 %>% filter(func_group == "Planktivore")

# pisc <- uvc.fish.2 %>% filter(genus_species == "cephalopholis argus")
# inve <- uvc.fish.2 %>% filter(genus_species == "parupeneus multifasciatus")
# herb <- uvc.fish.2 %>% filter(genus_species == "scarus niger")
# detr <- uvc.fish.2 %>% filter(genus_species == "ctenochaetus striatus")
# plan <- uvc.fish.2 %>% filter(genus_species == "naso vlamingii")

# create a list of parameters for input into the negll.PLB and pPLB functions
set.params = function(df){
  # Creates a list of parameters for input into the negll.PLB and pPLB functions
  size <- df$size_cm
  log.size <- log(size)
  sum.log.size <- sum(log.size)
  min.size <- min(df$size_cm)
  max.size <- max(df$size_cm)
  out.list <- (list(size, log.size, sum.log.size, min.size, max.size))
  names(out.list) <- c("size", "log.size", "sum.log.size", "min.size", "max.size")
  return(out.list)
}

# inputs for each functional group
pisc.input <- set.params(pisc)
inve.input <- set.params(inve)
herb.input <- set.params(herb)
detr.input <- set.params(detr)
plan.input <- set.params(plan)

# MLE settings
mgpVals <- c(1.6,0.5,0) # mgp values   2.0, 0.5, 0
xLim <- 10^par("usr")[1:2]
yLim <- 10^par("usr")[3:4]

# MLE size spectrum functions ---------------------------------------------
mle_b = function(x, log_x, sum_log_x, x_min, x_max){ # function added by PC
    # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
    # as a starting point for nlm for MLE of b for PLB model.
    PL.bMLE = 1/( log(min(x)) - sum_log_x/length(x)) - 1
    PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
        xmin=x_min, xmax=x_max, sumlogx=sum_log_x) #, print.level=2 )
    PLB.bMLE = PLB.minLL$estimate
    PLB.return = list(PLB.bMLE, PLB.minLL)
    return(PLB.return)
}

negLL.PLB = function(b, x, n, xmin, xmax, sumlogx){
  # Calculates the negative log-likelihood of the parameters b, xmin and xmax
  #  given data x for the PLB model. Returns the negative log-likelihood. Will
  #  be called by nlm or similar, but xmin and xmax are just estimated as the
  #  min and max of the data, not numerically using likelihood.
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood
  #   x: vector of values of data (e.g. masses of individual fish)
  #   n: length(x), have as an input to avoid repeatedly calculating it
  #   xmin: minimum value of x, have as an input to avoid repeatedly calculating
  #   xmax: maximum value of x, have as an input to avoid repeatedly calculating
  #   sumlogx: sum(log(x)) as an input, to avoid repeatedly calculating
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in negLL.PLB")
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumlogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumlogx
      }
    return(neglogLL)
}

pPLB = function(x = 10, b = -2, xmin = 1, xmax = 100){
  # Computes probability distribution function, P(X <= x),  for a
  #   bounded power-law (Pareto) distribution
  #
  # Args:
  #   x: vector of values at which to compute the distribution function
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  # Returns:
  #   vector of probability distribution values P(X <= x) corresponding to x
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y = 0 * x     # so have zeros where x < xmin
    y[x > xmax] = 1  # 1 for x > xmax
    if(b != -1)
        {  xmintobplus1 = xmin^(b+1)
           denom = xmax^(b+1) - xmintobplus1
           y[x >= xmin & x <= xmax] =
               ( x[x >= xmin & x <= xmax]^(b + 1) - xmintobplus1 ) / denom
        } else
        {  logxmin = log(xmin)
           denom = log(xmax) - logxmin
           y[x >= xmin & x <= xmax] =
               ( log( x[x >= xmin & x <= xmax] ) - logxmin ) / denom
        }
    return(y)
  }

# MLE Piscivores ----------------------------------------------------------
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.pisc <- mle_b(x=pisc.input$size, log_x=pisc.input$log.size, sum_log_x=pisc.input$sum.log.size,
                         x_min=pisc.input$min.size, x_max=pisc.input$max.size)
PLB.bMLE.pisc.b <- PLB.return.pisc[[1]]
PLB.minLL.pisc.b <- PLB.return.pisc[[2]]
# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.pisc.b <- PLB.minLL.pisc.b$minimum
x <- pisc.input$size
# plot(sort(pisc.input$size, decreasing=TRUE), 1:length(pisc.input$size), log="xy",
#      ylab = expression(paste("Number of body sizes", " ">=" ", italic("x"))), mgp=mgpVals,
#      xlab=expression(paste("Body sizes, ", italic(x), " (cm)")),
#      xlim = c(pisc.input$min.size, pisc.input$max.size), ylim = c(1, length(pisc.input$size)), axes=FALSE)
#      logTicks(xLim, yLim, xLabelBig = c(0, 1, 5, 10))   # Tick marks
#      x.PLB = seq(min(pisc.input$size), max(pisc.input$size), length=1000) # x values to plot PLB. Note
#                                                                           # that these encompass the data, and are not based
#                                                                           # on the binning (in MEE Figure 6 the line starts as
#                                                                           # min(x), not the first bin.
#      y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.pisc.b, xmin = min(x.PLB), xmax = max(x.PLB))) * length(pisc.input$size)
#      lines(x.PLB, y.PLB, col = "red", lwd = 2)
#      text(x = 16.5, y = 5, labels = "Piscivores", cex = 1.1, pos = 1, col = "black")
#      spectra.text <- as.character(round(PLB.bMLE.pisc.b, 2))
#      text(x = 16, y = 2.5, labels = bquote(paste(italic("b = "), .(spectra.text))), cex = 1.1, pos = 1, col="black")
# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.pisc.b - 0.5, PLB.bMLE.pisc.b + 0.5, 0.00001)
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
     PLB.LLvals[i] = negLL.PLB(bvec[i], x=pisc.input$size, n=length(pisc.input$size), xmin=pisc.input$min.size,
     xmax=pisc.input$max.size, sumlogx=pisc.input$sum.log.size)
}
critVal = PLB.minNegLL.pisc.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
PLB.bMLE.pisc.CI <- c(min(bIn95), max(bIn95))
# # To add just the curves at the limits of the 95% confidence interval of b:
# for(i in c(1, length(bIn95))){
#      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
#      xmax = max(x.PLB))) * length(pisc.input$size), col="red", lty=2)
# }
pisc.out <- list(PLB.bMLE.pisc.b, PLB.bMLE.pisc.CI)

# MLE Invertivores ----------------------------------------------------------
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.inve <- mle_b(x=inve.input$size, log_x=inve.input$log.size, sum_log_x=inve.input$sum.log.size,
                         x_min=inve.input$min.size, x_max=inve.input$max.size)
PLB.bMLE.inve.b <- PLB.return.inve[[1]]
PLB.minLL.inve.b <- PLB.return.inve[[2]]
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.inve.b <- PLB.minLL.inve.b$minimum
# x <- inve.input$size
# plot(sort(inve.input$size, decreasing=TRUE), 1:length(inve.input$size), log="xy",
#      ylab = expression(paste("Number of body sizes", " ">=" ", italic("x"))), mgp=mgpVals,
#      xlab=expression(paste("Body sizes, ", italic(x), " (cm)")),
#      xlim = c(inve.input$min.size, inve.input$max.size), ylim = c(1, length(inve.input$size)), axes=FALSE)
#      logTicks(xLim, yLim, xLabelBig = c(0, 1, 5, 10))   # Tick marks
#      x.PLB = seq(min(inve.input$size), max(inve.input$size), length=1000) # x values to plot PLB. Note
#                                                                           # that these encompass the data, and are not based
#                                                                           # on the binning (in MEE Figure 6 the line starts as
#                                                                           # min(x), not the first bin.
#      y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.inve.b, xmin = min(x.PLB), xmax = max(x.PLB))) * length(inve.input$size)
#      lines(x.PLB, y.PLB, col = "red", lwd = 2)
#      text(x = 16.5, y = 5, labels = "Invertivores", cex = 1.1, pos = 1, col = "black")
#      spectra.text <- as.character(round(PLB.bMLE.inve.b, 2))
#      text(x = 16, y = 2.5, labels = bquote(paste(italic("b = "), .(spectra.text))), cex = 1.1, pos = 1, col="black")

# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.inve.b - 0.5, PLB.bMLE.inve.b + 0.5, 0.00001)
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
     PLB.LLvals[i] = negLL.PLB(bvec[i], x=inve.input$size, n=length(inve.input$size), xmin=inve.input$min.size,
     xmax=inve.input$max.size, sumlogx=inve.input$sum.log.size)
}
critVal = PLB.minNegLL.inve.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
PLB.bMLE.inve.CI <- c(min(bIn95), max(bIn95))
# # To add just the curves at the limits of the 95% confidence interval of b:
# for(i in c(1, length(bIn95))){
#      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
#      xmax = max(x.PLB))) * length(inve.input$size), col="red", lty=2)
# }
inve.out <- list(PLB.bMLE.inve.b, PLB.bMLE.inve.CI)

# MLE Herbivores ----------------------------------------------------------
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.herb <- mle_b(x=herb.input$size, log_x=herb.input$log.size, sum_log_x=herb.input$sum.log.size,
                         x_min=herb.input$min.size, x_max=herb.input$max.size)
PLB.bMLE.herb.b <- PLB.return.herb[[1]]
PLB.minLL.herb.b <- PLB.return.herb[[2]]
# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.herb.b <- PLB.minLL.herb.b$minimum
x <- herb.input$size
# plot(sort(herb.input$size, decreasing=TRUE), 1:length(herb.input$size), log="xy",
#      ylab = expression(paste("Number of body sizes", " ">=" ", italic("x"))), mgp=mgpVals,
#      xlab=expression(paste("Body sizes, ", italic(x), " (cm)")),
#      xlim = c(herb.input$min.size, herb.input$max.size), ylim = c(1, length(herb.input$size)), axes=FALSE)
#      logTicks(xLim, yLim, xLabelBig = c(0, 1, 5, 10))   # Tick marks
#      x.PLB = seq(min(herb.input$size), max(herb.input$size), length=1000) # x values to plot PLB. Note
#                                                                           # that these encompass the data, and are not based
#                                                                           # on the binning (in MEE Figure 6 the line starts as
#                                                                           # min(x), not the first bin.
#      y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.herb.b, xmin = min(x.PLB), xmax = max(x.PLB))) * length(herb.input$size)
#      lines(x.PLB, y.PLB, col = "red", lwd = 2)
#      text(x = 16.5, y = 5, labels = "Herbivores", cex = 1.1, pos = 1, col = "black")
#      spectra.text <- as.character(round(PLB.bMLE.herb.b, 2))
#      text(x = 16, y = 2.5, labels = bquote(paste(italic("b = "), .(spectra.text))), cex = 1.1, pos = 1, col="black")

# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.herb.b - 0.5, PLB.bMLE.herb.b + 0.5, 0.00001)
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
     PLB.LLvals[i] = negLL.PLB(bvec[i], x=herb.input$size, n=length(herb.input$size), xmin=herb.input$min.size,
     xmax=herb.input$max.size, sumlogx=herb.input$sum.log.size)
}
critVal = PLB.minNegLL.herb.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
PLB.bMLE.herb.CI <- c(min(bIn95), max(bIn95))
# # To add just the curves at the limits of the 95% confidence interval of b:
# for(i in c(1, length(bIn95))){
#      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
#      xmax = max(x.PLB))) * length(herb.input$size), col="red", lty=2)
# }
herb.out <- list(PLB.bMLE.herb.b, PLB.bMLE.herb.CI)

# MLE Detritivores ----------------------------------------------------------
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.detr <- mle_b(x=detr.input$size, log_x=detr.input$log.size, sum_log_x=detr.input$sum.log.size,
                         x_min=detr.input$min.size, x_max=detr.input$max.size)
PLB.bMLE.detr.b <- PLB.return.detr[[1]]
PLB.minLL.detr.b <- PLB.return.detr[[2]]
# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.detr.b <- PLB.minLL.detr.b$minimum
x <- detr.input$size
# plot(sort(detr.input$size, decreasing=TRUE), 1:length(detr.input$size), log="xy",
#      ylab = expression(paste("Number of body sizes", " ">=" ", italic("x"))), mgp=mgpVals,
#      xlab=expression(paste("Body sizes, ", italic(x), " (cm)")),
#      xlim = c(detr.input$min.size, detr.input$max.size), ylim = c(1, length(detr.input$size)), axes=FALSE)
#      logTicks(xLim, yLim, xLabelBig = c(0, 1, 5, 10))   # Tick marks
#      x.PLB = seq(min(detr.input$size), max(detr.input$size), length=1000) # x values to plot PLB. Note
#                                                                           # that these encompass the data, and are not based
#                                                                           # on the binning (in MEE Figure 6 the line starts as
#                                                                           # min(x), not the first bin.
#      y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.detr.b, xmin = min(x.PLB), xmax = max(x.PLB))) * length(detr.input$size)
#      lines(x.PLB, y.PLB, col = "red", lwd = 2)
#      text(x = 13, y = 5, labels = "Detritivores", cex = 1.1, pos = 1, col = "black")
#      spectra.text <- as.character(round(PLB.bMLE.detr.b, 2))
#      text(x = 12.8, y = 2.5, labels = bquote(paste(italic("b = "), .(spectra.text))), cex = 1.1, pos = 1, col="black")

# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.detr.b - 0.5, PLB.bMLE.detr.b + 0.5, 0.00001)
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
     PLB.LLvals[i] = negLL.PLB(bvec[i], x=detr.input$size, n=length(detr.input$size), xmin=detr.input$min.size,
     xmax=detr.input$max.size, sumlogx=detr.input$sum.log.size)
}
critVal = PLB.minNegLL.detr.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
PLB.bMLE.detr.CI <- c(min(bIn95), max(bIn95))
# # To add just the curves at the limits of the 95% confidence interval of b:
# for(i in c(1, length(bIn95))){
#      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
#      xmax = max(x.PLB))) * length(detr.input$size), col="red", lty=2)
# }
detr.out <- list(PLB.bMLE.detr.b, PLB.bMLE.detr.CI)

# MLE Planktivores ----------------------------------------------------------
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.plan <- mle_b(x=plan.input$size, log_x=plan.input$log.size, sum_log_x=plan.input$sum.log.size,
                         x_min=plan.input$min.size, x_max=plan.input$max.size)
PLB.bMLE.plan.b <- PLB.return.plan[[1]]
PLB.minLL.plan.b <- PLB.return.plan[[2]]
# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.plan.b <- PLB.minLL.plan.b$minimum
x <- plan.input$size
# plot(sort(plan.input$size, decreasing=TRUE), 1:length(plan.input$size), log="xy",
#      ylab = expression(paste("Number of body sizes", " ">=" ", italic("x"))), mgp=mgpVals,
#      xlab=expression(paste("Body sizes, ", italic(x), " (cm)")),
#      xlim = c(plan.input$min.size, plan.input$max.size), ylim = c(1, length(plan.input$size)), axes=FALSE)
#      logTicks(xLim, yLim, xLabelBig = c(0, 1, 5, 10))   # Tick marks
#      x.PLB = seq(min(plan.input$size), max(plan.input$size), length=1000) # x values to plot PLB. Note
#                                                                           # that these encompass the data, and are not based
#                                                                           # on the binning (in MEE Figure 6 the line starts as
#                                                                           # min(x), not the first bin.
#      y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.plan.b, xmin = min(x.PLB), xmax = max(x.PLB))) * length(plan.input$size)
#      lines(x.PLB, y.PLB, col = "red", lwd = 2)
#      text(x = 13, y = 5, labels = "Planktivores", cex = 1.1, pos = 1, col = "black")
#      spectra.text <- as.character(round(PLB.bMLE.plan.b, 2))
#      text(x = 12.8, y = 2.5, labels = bquote(paste(italic("b = "), .(spectra.text))), cex = 1.1, pos = 1, col="black")

# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.plan.b - 0.5, PLB.bMLE.plan.b + 0.5, 0.00001)
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
     PLB.LLvals[i] = negLL.PLB(bvec[i], x=plan.input$size, n=length(plan.input$size), xmin=plan.input$min.size,
     xmax=plan.input$max.size, sumlogx=plan.input$sum.log.size)
}
critVal = PLB.minNegLL.plan.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
PLB.bMLE.plan.CI <- c(min(bIn95), max(bIn95))
# # To add just the curves at the limits of the 95% confidence interval of b:
# for(i in c(1, length(bIn95))){
#      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
#      xmax = max(x.PLB))) * length(plan.input$size), col="red", lty=2)
# }
plan.out <- list(PLB.bMLE.plan.b, PLB.bMLE.plan.CI)
