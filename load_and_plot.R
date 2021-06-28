# ------------------------------------------------------------------------------------------------------------------
## File 'load_and_plot.R'
## Gear-based fisheries population model
##
## Authors: Paul Carvalho
##          Austin Humphries
##
## Description: Create manuscript plots. Coral reef fisheries model of nine functional groups to test biomass
##              and catch response to various gear-based management scenarios. The model losely represents 
##			 the Wakatobi coral reef fishery in Indonesia, where certain model parameters were derived.

# Clean workspace
rm(list = ls())

# --------------------------------------------------- DIRECTORIES ---------------------------------------------------

# setwd("C:/Users/pgcar/Google Drive/Paul Carvalho/dissertation/chapter 4/model") # PC
# setwd("~/Google Drive/grad students/Paul Carvalho/dissertation/chapter 4/model")  # AH

# --------------------------------------------------- LIBRARIES ---------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggpubr)
library(pracma)

# --------------------------------------------------- LOAD MODEL DATA ---------------------------------------------------

load("model_data.RData")

# --------------------------------------------------- FIG 1: total biomass and catch ---------------------------------------------------

# Colorblind friendly palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# Create dataframe for total biomass and catch of all scenarios
scenarios <- c("All gears", "Hook-and-line ban", "Net ban", "Spear ban", "Hook-and-line", "Net", "Spear")
plot1.df  <- data.frame(scenario = rep(scenarios, each = length(effort)),
				    effort   = rep(effort, length(scenarios)),
				    biomass  = c(scen.1.out[[1]], scen.2.out[[1]], scen.3.out[[1]], scen.4.out[[1]], scen.5.out[[1]], scen.6.out[[1]], scen.7.out[[1]]),
				    catch    = c(scen.1.out[[2]], scen.2.out[[2]], scen.3.out[[2]], scen.4.out[[2]], scen.5.out[[2]], scen.6.out[[2]], scen.7.out[[2]]))
bs.index <- c(6, 10, 15)
plot1.ci <- data.frame(scenario = rep(scenarios, each = 3),
                       effort   = rep(effort.1.bs, length(scenarios)),
                       B.mu     = c(scen.1.out[[1]][bs.index], scen.2.out[[1]][bs.index], scen.3.out[[1]][bs.index], scen.4.out[[1]][bs.index], scen.5.out[[1]][bs.index], scen.6.out[[1]][bs.index], scen.7.out[[1]][bs.index]),
                       CB.mu    = c(scen.1.out[[2]][bs.index], scen.2.out[[2]][bs.index], scen.3.out[[2]][bs.index], scen.4.out[[2]][bs.index], scen.5.out[[2]][bs.index], scen.6.out[[2]][bs.index], scen.7.out[[2]][bs.index]),
                       B.lo     = c(scen.1.bs[[1]]$B.lo - (scen.1.bs[[1]]$B.mu - scen.1.out[[1]][bs.index]), scen.2.bs[[1]]$B.lo - (scen.2.bs[[1]]$B.mu - scen.2.out[[1]][bs.index]), scen.3.bs[[1]]$B.lo - (scen.3.bs[[1]]$B.mu - scen.3.out[[1]][bs.index]), scen.4.bs[[1]]$B.lo - (scen.4.bs[[1]]$B.mu - scen.4.out[[1]][bs.index]), scen.5.bs[[1]]$B.lo - (scen.5.bs[[1]]$B.mu - scen.5.out[[1]][bs.index]), scen.6.bs[[1]]$B.lo - (scen.6.bs[[1]]$B.mu - scen.6.out[[1]][bs.index]), scen.7.bs[[1]]$B.lo - (scen.7.bs[[1]]$B.mu - scen.7.out[[1]][bs.index])),
                       B.up     = c(scen.1.bs[[1]]$B.up - (scen.1.bs[[1]]$B.mu - scen.1.out[[1]][bs.index]), scen.2.bs[[1]]$B.up - (scen.2.bs[[1]]$B.mu - scen.2.out[[1]][bs.index]), scen.3.bs[[1]]$B.up - (scen.3.bs[[1]]$B.mu - scen.3.out[[1]][bs.index]), scen.4.bs[[1]]$B.up - (scen.4.bs[[1]]$B.mu - scen.4.out[[1]][bs.index]), scen.5.bs[[1]]$B.up - (scen.5.bs[[1]]$B.mu - scen.5.out[[1]][bs.index]), scen.6.bs[[1]]$B.up - (scen.6.bs[[1]]$B.mu - scen.6.out[[1]][bs.index]), scen.7.bs[[1]]$B.up - (scen.7.bs[[1]]$B.mu - scen.7.out[[1]][bs.index])),                       
                       CB.lo    = c(scen.1.bs[[1]]$CB.lo - (scen.1.bs[[1]]$CB.mu - scen.1.out[[2]][bs.index]), scen.2.bs[[1]]$CB.lo - (scen.2.bs[[1]]$CB.mu - scen.2.out[[2]][bs.index]), scen.3.bs[[1]]$CB.lo - (scen.3.bs[[1]]$CB.mu - scen.3.out[[2]][bs.index]), scen.4.bs[[1]]$CB.lo - (scen.4.bs[[1]]$CB.mu - scen.4.out[[2]][bs.index]), scen.5.bs[[1]]$CB.lo - (scen.5.bs[[1]]$CB.mu - scen.5.out[[2]][bs.index]), scen.6.bs[[1]]$CB.lo - (scen.6.bs[[1]]$CB.mu - scen.6.out[[2]][bs.index]), scen.7.bs[[1]]$CB.lo - (scen.7.bs[[1]]$CB.mu - scen.7.out[[2]][bs.index])),
                       CB.up    = c(scen.1.bs[[1]]$CB.up - (scen.1.bs[[1]]$CB.mu - scen.1.out[[2]][bs.index]), scen.2.bs[[1]]$CB.up - (scen.2.bs[[1]]$CB.mu - scen.2.out[[2]][bs.index]), scen.3.bs[[1]]$CB.up - (scen.3.bs[[1]]$CB.mu - scen.3.out[[2]][bs.index]), scen.4.bs[[1]]$CB.up - (scen.4.bs[[1]]$CB.mu - scen.4.out[[2]][bs.index]), scen.5.bs[[1]]$CB.up - (scen.5.bs[[1]]$CB.mu - scen.5.out[[2]][bs.index]), scen.6.bs[[1]]$CB.up - (scen.6.bs[[1]]$CB.mu - scen.6.out[[2]][bs.index]), scen.7.bs[[1]]$CB.up - (scen.7.bs[[1]]$CB.mu - scen.7.out[[2]][bs.index])))

# Convert to relative values
plot1.df <- plot1.df %>%
	dplyr::mutate(effort.rel = effort/max(effort)) %>%    # relative effort
     dplyr::mutate(biomass.rel = biomass/max(biomass)) %>% # relative biomass
     dplyr::mutate(catch.rel = catch/max(catch)) %>%       # relative catch
     dplyr::mutate(sq.diff = (biomass.rel-0.5)^2)          # find squared difference between relative biomass and 0.5 to identify effort and catch at 0.5B0
plot1.ci <- plot1.ci %>%
     dplyr::mutate(effort.rel = effort/max(plot1.df$effort)) %>% # relative effort
     dplyr::mutate(B.mu = B.mu/max(plot1.df$biomass)) %>%
     dplyr::mutate(B.lo = B.lo/max(plot1.df$biomass)) %>%
     dplyr::mutate(B.up = B.up/max(plot1.df$biomass)) %>%
     dplyr::mutate(CB.mu = CB.mu/max(plot1.df$catch)) %>%
     dplyr::mutate(CB.lo = CB.lo/max(plot1.df$catch)) %>%
     dplyr::mutate(CB.up = CB.up/max(plot1.df$catch))

# Find catch and effort at 0.5B0
plot1Half.df <- plot1.df %>%
	dplyr::group_by(scenario) %>%
	dplyr::summarise(min = min(sq.diff)) %>%
	left_join(., plot1.df, by = c("min" = "sq.diff")) %>%
	dplyr::select(scenario = scenario.x, effort.rel, catch.rel)

# Plot
plot1.B <- ggplot() +
	geom_line(data = plot1.df, aes(x = effort.rel, y = biomass.rel, color = scenario), lwd = 0.5) +
	geom_hline(yintercept = 0.5, lty = "dashed", lwd = 1, alpha = 0.3) +
     geom_point(aes(x = plot1.ci$effort.rel, y = plot1.ci$B.mu, color = plot1.ci$scenario), size = 0.75) +     
     geom_errorbar(aes(x = plot1.ci$effort.rel, ymin = plot1.ci$B.lo, ymax = plot1.ci$B.up, color = plot1.ci$scenario), lwd = 0.5, width = 0.01) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
	labs(x = "", y = expression(paste("Stock biomass (relative to ", B[paste(" ", italic("F"),"=0")], ")"))) +
	theme_classic() +
	theme(legend.title = element_blank(), legend.position = c(0.85,0.8), plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))

plot1.C <- ggplot() +
	geom_line(data = plot1.df, aes(x = effort.rel, y = catch.rel, color = scenario), lwd = 0.5) +
     geom_point(aes(x = plot1.ci$effort.rel, y = plot1.ci$CB.mu, color = plot1.ci$scenario), size = 0.75) +
     geom_errorbar(aes(x = plot1.ci$effort.rel, ymin = plot1.ci$CB.lo, ymax = plot1.ci$CB.up, color = plot1.ci$scenario), lwd = 0.5, width = 0.01) +
     scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
	labs(x = "Effort", y = expression(paste("Catch (relative to ", C[italic("max")], " across all scenarios)"))) +
	theme_classic() +
	theme(legend.position = "none", plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
	
ggarrange(plot1.B, plot1.C, nrow=2, ncol=1, labels=c("A","B"))


# --------------------------------------------------- FIG 2: mean length for stock and catch --------------------------------------------------

# Note: Relative mean weight and length for biomass and catch were equivalent and, thus, we only included mean length.
# Create dataframe for mean length (biomass and catch)
plot2.df <- data.frame(scenario = rep(scenarios, each = length(effort)),
				   effort = rep(effort, length(scenarios)),
				   BmeanL = c(scen.1.out[[6]], scen.2.out[[6]], scen.3.out[[6]], scen.4.out[[6]], scen.5.out[[6]], scen.6.out[[6]], scen.7.out[[6]]),
				   CmeanL = c(scen.1.out[[9]], scen.2.out[[9]], scen.3.out[[9]], scen.4.out[[9]], scen.5.out[[9]], scen.6.out[[9]], scen.7.out[[9]]))

# Convert to relative values
plot2.df <- plot2.df %>%
	mutate(effort.rel = effort/max(effort)) %>%               # relative effort
	mutate(BmeanL.rel = BmeanL/max(BmeanL)) %>%			   # relative mean length (biomass)
	mutate(CmeanL.rel = CmeanL/max(CmeanL, na.rm = TRUE))     # relative mean length (catch)
plot2.half <- plot1Half.df %>%
	left_join(., plot2.df, by = c("effort.rel", "scenario"))

# Plot
plot2.BmeanL <- ggplot() +
	geom_line(data = plot2.df, aes(x = effort.rel, y = BmeanL.rel, color = scenario), lwd = 1) +
	geom_point(data = plot2.half, aes(x = effort.rel, y = BmeanL.rel, color = scenario), size = 2.5, alpha = 0.75) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
	labs(x = "", y = "Relative mean length (stock)") +
	theme_classic() +
	theme(legend.title = element_blank(), legend.position = c(0.85,0.8), plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
plot2.CmeanL <- ggplot() +
	geom_line(data = plot2.df, aes(x = effort.rel, y = CmeanL.rel, color = scenario), lwd = 1) +
	geom_point(data = plot2.half, aes(x = effort.rel, y = CmeanL.rel, color = scenario), size = 2.5, alpha = 0.75) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
	labs(x = "Effort", y = "Relative mean length (catch)") +
	theme_classic() +
	theme(legend.title = element_blank(), legend.position = "none", plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))

ggarrange(plot2.BmeanL, plot2.CmeanL, nrow=2, ncol=1, labels=c("A","B"))

# --------------------------------------------------- FIG 3-4: functional group biomass and catch ---------------------------------------------------

# Create dataframe with functional group biomass and catch for all scenarios
fg <- c("Browsers", "Detritivores", "Excavators/scrapers", "Grazers", "Macro-invertivores", "Micro-invertivores", "Pisci-invertivores", "Piscivores", "Planktivores") 
plot3.df <- data.frame(scenario = rep(scenarios, each = (length(fg) * length(effort))),
				   fg = rep(fg, each = length(effort), times = length(scenarios)),
				   effort = rep(effort, times = (length(fg) * length(scenarios))),
				   biomass = c(as.vector(t(scen.1.out[[3]])), as.vector(t(scen.2.out[[3]])), as.vector(t(scen.3.out[[3]])), as.vector(t(scen.4.out[[3]])), as.vector(t(scen.5.out[[3]])), as.vector(t(scen.6.out[[3]])), as.vector(t(scen.7.out[[3]]))),
				   catch = c(as.vector(t(scen.1.out[[4]])), as.vector(t(scen.2.out[[4]])), as.vector(t(scen.3.out[[4]])), as.vector(t(scen.4.out[[4]])), as.vector(t(scen.5.out[[4]])), as.vector(t(scen.6.out[[4]])), as.vector(t(scen.7.out[[4]]))))
plot3.ci <- data.frame(scenario = rep(scenarios, each = (length(fg) * length(bs.index))),
                       fg = rep(fg, each = length(bs.index), times = length(scenarios)),
                       effort = rep(effort[bs.index], times = (length(fg) * length(scenarios))),
                       biomass = c(as.vector(t(scen.1.out[[3]][,bs.index])), as.vector(t(scen.2.out[[3]][,bs.index])), as.vector(t(scen.3.out[[3]][,bs.index])), as.vector(t(scen.4.out[[3]][,bs.index])), as.vector(t(scen.5.out[[3]][,bs.index])), as.vector(t(scen.6.out[[3]][,bs.index])), as.vector(t(scen.7.out[[3]][,bs.index]))),
                       catch = c(as.vector(t(scen.1.out[[4]][,bs.index])), as.vector(t(scen.2.out[[4]][,bs.index])), as.vector(t(scen.3.out[[4]][,bs.index])), as.vector(t(scen.4.out[[4]][,bs.index])), as.vector(t(scen.5.out[[4]][,bs.index])), as.vector(t(scen.6.out[[4]][,bs.index])), as.vector(t(scen.7.out[[4]][,bs.index]))),
                       B.up = c(scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$B.up - (scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$B.mu - as.vector(t(scen.1.out[[3]][,bs.index]))), scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$B.up - (scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$B.mu - as.vector(t(scen.2.out[[3]][,bs.index]))), scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$B.up - (scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$B.mu - as.vector(t(scen.3.out[[3]][,bs.index]))), scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$B.up - (scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$B.mu - as.vector(t(scen.4.out[[3]][,bs.index]))), scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$B.up - (scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$B.mu - as.vector(t(scen.5.out[[3]][,bs.index]))), scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$B.up - (scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$B.mu - as.vector(t(scen.6.out[[3]][,bs.index]))), scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$B.up - (scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$B.mu - as.vector(t(scen.7.out[[3]][,bs.index])))),
                       B.lo = c(scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$B.lo - (scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$B.mu - as.vector(t(scen.1.out[[3]][,bs.index]))), scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$B.lo - (scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$B.mu - as.vector(t(scen.2.out[[3]][,bs.index]))), scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$B.lo - (scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$B.mu - as.vector(t(scen.3.out[[3]][,bs.index]))), scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$B.lo - (scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$B.mu - as.vector(t(scen.4.out[[3]][,bs.index]))), scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$B.lo - (scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$B.mu - as.vector(t(scen.5.out[[3]][,bs.index]))), scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$B.lo - (scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$B.mu - as.vector(t(scen.6.out[[3]][,bs.index]))), scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$B.lo - (scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$B.mu - as.vector(t(scen.7.out[[3]][,bs.index])))),
                       CB.up = c(scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$CB.up - (scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.1.out[[4]][,bs.index]))), scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$CB.up - (scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.2.out[[4]][,bs.index]))), scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$CB.up - (scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.3.out[[4]][,bs.index]))), scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$CB.up - (scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.4.out[[4]][,bs.index]))), scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$CB.up - (scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.5.out[[4]][,bs.index]))), scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$CB.up - (scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.6.out[[4]][,bs.index]))), scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$CB.up - (scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.7.out[[4]][,bs.index])))),
                       CB.lo = c(scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$CB.lo - (scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.1.out[[4]][,bs.index]))), scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$CB.lo - (scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.2.out[[4]][,bs.index]))), scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$CB.lo - (scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.3.out[[4]][,bs.index]))), scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$CB.lo - (scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.4.out[[4]][,bs.index]))), scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$CB.lo - (scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.5.out[[4]][,bs.index]))), scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$CB.lo - (scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.6.out[[4]][,bs.index]))), scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$CB.lo - (scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$CB.mu - as.vector(t(scen.7.out[[4]][,bs.index])))))

# force negative confidence intervals to 0
plot3.ci$B.lo[which(plot3.ci$B.lo < 0)]   <- 0
plot3.ci$CB.lo[which(plot3.ci$CB.lo < 0)] <- 0

# Iterate through functional groups to construct plots
plot3.B <- vector(mode="list", length=length(fg)) # empty list to store ggplots
plot4.C <- vector(mode="list", length=length(fg)) # empty list to store ggplots

for(i in 1:length(fg)){
   fg.i <- fg[i] # Save name of functional group
   
   # Create dataframe for fg.i
   fg.df <- plot3.df %>% 
   	filter(fg == fg.i) %>%                         # get data for fg.i
   	mutate(effort.rel = effort/max(effort)) %>%    # relative effort
	mutate(biomass.rel = biomass/max(biomass)) %>% # relative biomass
   	mutate(catch.rel = catch/max(catch))           # relative catch
   
   # Create dataframe for confidence intervals
   fg.ci <- plot3.ci %>%
     filter(fg == fg.i) %>%
     mutate(effort.rel = effort/max(fg.df$effort)) %>%
     mutate(biomass.rel = biomass/max(fg.df$biomass)) %>%
     mutate(catch.rel = catch/max(fg.df$catch)) %>%
     mutate(B.up = B.up/max(fg.df$biomass)) %>%
     mutate(B.lo = B.lo/max(fg.df$biomass)) %>%
     mutate(CB.up = CB.up/max(fg.df$catch)) %>%
     mutate(CB.lo = CB.lo/max(fg.df$catch))
   
   tmp.df <- plot1Half.df %>% 
   	left_join(., fg.df, by = c("effort.rel","scenario")) # get biomass and catch of fg.1 when total biomass = 0.5B0
   
   # create biomass plot for fg.i and all scenarios
   fg.plotB <- ggplot() +
   	geom_line(data = fg.df, aes(x = effort.rel, y = biomass.rel, color = scenario), lwd = 0.35) +
     geom_point(data = fg.ci, aes(x = effort.rel, y = biomass.rel, color = scenario), size = 0.25) +
     geom_errorbar(data = fg.ci, aes(x = effort.rel, ymin = B.lo, ymax = B.up, color = scenario), lwd = 0.35, width = 0.02) +
   	scale_y_continuous(expand = c(0,0), limits = c(0, 1.02)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
   	labs(x = "", y = "", title = fg.i) +
   	theme_classic() +
   	theme(plot.title = element_text(hjust = 0.5, size = 9), legend.title = element_blank(), legend.position = "none", plot.margin = margin(0.4,0.4,0.4,0.4,"cm"))
   	 
   plot3.B[[i]] <- fg.plotB
   
   # create catch plot for fg.i and all scenarios
   fg.plotC <- ggplot() +
   	geom_line(data = fg.df, aes(x = effort.rel, y = catch.rel, color = scenario), lwd = 0.35) +
     geom_point(data = fg.ci, aes(x = effort.rel, y = catch.rel, color = scenario), size = 0.25) +
     geom_errorbar(data = fg.ci, aes(x = effort.rel, ymin = CB.lo, ymax = CB.up, color = scenario), lwd = 0.35, width = 0.02) +
   	scale_y_continuous(expand = c(0,0), limits = c(0,1.02)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
   	labs(x = "", y = "", title = fg.i) +
   	theme_classic() +
   	theme(plot.title = element_text(hjust = 0.5, size = 9), legend.title = element_blank(), legend.position = "none", plot.margin = margin(0.4,0.4,0.4,0.4,"cm"))
   
   plot4.C[[i]] <- fg.plotC
}

plot3 <- ggarrange(plot3.B[[1]],plot3.B[[2]],plot3.B[[3]],plot3.B[[4]],plot3.B[[5]],
			    plot3.B[[6]],plot3.B[[7]],plot3.B[[8]],plot3.B[[9]], nrow=3, ncol=3,
			    labels=c("A","B","C","D","E","F","G","H","I"), common.legend = TRUE, legend = "bottom",
			    label.x = 0.1)

annotate_figure(plot3,
			 bottom = text_grob("Effort", vjust=-6),
			 left = text_grob(expression(paste("Stock biomass (relative to ", B[paste(" ", italic("F"),"=0")], ")")), rot = 90, vjust = 1.5, hjust = 0.15))

plot4 <- ggarrange(plot4.C[[1]],plot4.C[[2]],plot4.C[[3]],plot4.C[[4]],plot4.C[[5]],
			    plot4.C[[6]],plot4.C[[7]],plot4.C[[8]],plot4.C[[9]], nrow=3, ncol=3,
		         labels=c("A","B","C","D","E","F","G","H","I"), common.legend = TRUE, legend = "bottom",
			    label.x = 0.1)

annotate_figure(plot4,
			 bottom = text_grob("Effort", vjust=-6),
			 left = text_grob(expression(paste("Catch (relative to ", C[italic(" max")], ")")), rot = 90, vjust = 1.5, hjust = 0.15))

# --------------------------------------------------- FIG SX: sensitivity of biomass and catch to removing functional groups ---------------------------------------------------

# Remove a functional group from total biomass
fg.less1 <- plot3.df %>%
	filter(fg != "Planktivores") %>%
	dplyr::group_by(scenario, effort) %>%
	dplyr::summarise(biomass = sum(biomass),
			       catch = sum(catch)) %>%
	ungroup() %>%
	mutate(effort.rel = effort/max(effort)) %>%
	mutate(bio.rel = biomass/max(biomass)) %>%
	mutate(catch.rel = catch/max(catch))

plot1.B

fg.less1.plot <- ggplot() +
	geom_line(data = fg.less1, aes(x = effort.rel, y = bio.rel, color = scenario), lwd = 1) +
	geom_hline(yintercept = 0.5, lty = "dashed", lwd = 1, alpha = 0.3) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
	labs(x = "", y = expression(paste("Stock biomass (relative to ", B[paste(" ", italic("F"),"=0")], ")"))) +
	theme_classic() +
	theme(legend.title = element_blank(),
		 legend.position = c(0.85,0.8),
		 plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))

ggarrange(plot1.B, fg.less1.plot)

# --------------------------------------------------- FIG SX: lognormal selectivity distributions ---------------------------------------------------
x <- seq(0,65,0.5)

line <- dlnorm(x, meanlog = 2.9140, sdlog = 0.4098)
net <- dlnorm(x, meanlog = 2.8821, sdlog = 0.3304)
spear <- dlnorm(x, meanlog = 3.0396, sdlog = 0.2646)

sel.df <- data.frame(size = rep(x,3),
				 gear = rep(c("Line","Net","Spear"), each = length(x)),
				 sel = c(line,net,spear))

sel.plot <- ggplot() +
	geom_line(data = sel.df, aes(x = size, y = sel, color = gear), lwd = 1) +
	theme_classic() +
	scale_color_manual(values = cbp2[1:3]) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	labs(x = "Size (cm)", y = "Selectivity") +
	theme(legend.title = element_blank(),
		 legend.position = c(0.75,0.75),
		 plot.margin = margin(0.3,0.3,0.3,0.3,"cm"))
	

# --------------------------------------------------- FIG SX: sensitivity to gear specification ---------------------------------------------------

# Create data frame for different gear specifications
gear.spec <- c("0-cm increase", "1-cm increase", "2-cm increase")
gear.spec.df <- data.frame(gear.spec = rep(gear.spec, each = length(effort)),
					  effort = rep(effort, length(gear.spec)),
					  biomass = c(scen.1.out[[1]], q1sa.out[[1]], q2sa.out[[1]]),
					  catch = c(scen.1.out[[2]], q1sa.out[[2]], q2sa.out[[2]]))

# Convert to relative values
gear.spec.df <- gear.spec.df %>%
	mutate(effort.rel = effort/max(effort)) %>% 	  # relative effort
	mutate(biomass.rel = biomass/max(biomass)) %>% # relative biomass
	mutate(catch.rel = catch/max(catch)) %>%	  # relative catch
	mutate(sq.diff = (biomass.rel - 0.5)^2)	       # find squared difference between relative biomass and 0.5 to identify effort and catch at 0.5B0
	
# Find catch and effort at 0.5B0
gear.spec.half.df <- gear.spec.df %>%
	group_by(gear.spec) %>%
	summarise(min = min(sq.diff)) %>%
	left_join(., gear.spec.df, by = c("min" = "sq.diff")) %>%
	dplyr::select(gear.spec = gear.spec.x, effort.rel, catch.rel)
	
# Plot
gearspec.B <- ggplot() +
	geom_line(data = gear.spec.df, aes(x = effort.rel, y = biomass.rel, lty = gear.spec), lwd = 1) +
	geom_hline(yintercept = 0.5, lty = "dashed", lwd = 1, alpha = 0.3) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	labs(x = "", y = expression(paste("Stock biomass (relative to ", B[paste(" ", italic("F"),"=0")], ")"))) +
	theme_classic() +
	theme(legend.title = element_blank(),
		 legend.position = c(0.85,0.8),
		 plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
gearspec.C <- ggplot() +
	geom_line(data = gear.spec.df, aes(x = effort.rel, y = catch.rel, lty = gear.spec), lwd = 1) +
	geom_point(data = gear.spec.half.df, aes(x = effort.rel, y = catch.rel), size = 2.5, alpha = 0.75) +
	geom_linerange(data = gear.spec.half.df, aes(x = effort.rel, ymax = catch.rel, ymin = 0, lty = gear.spec), lwd = 1, alpha = 0.5) +
	geom_linerange(data = gear.spec.half.df, aes(xmin = 0, xmax = effort.rel, y = catch.rel, lty = gear.spec), lwd = 1, alpha = 0.5) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	labs(x = "Effort", y = expression(paste("Catch (relative to ", C[italic("max")], " across all scenarios)"))) +
	theme_classic() +
	theme(legend.position = "none",
		 plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
	
ggarrange(gearspec.B, gearspec.C, nrow=2, ncol=1, labels=c("A","B"))


# --------------------------------------------------- FIG SX: heatmaps for catchability and selectivity ---------------------------------------------------

# create dataframe for ggplot heatmap (geom_tile)
sc <- paste("[", L.lower,"-",L.upper, ")", sep="")
heatmap.df <- data.frame(fg = rep(fg, each=length(sc)),
					sc = rep(sc, times = length(fg)),
					qs.hook = as.vector(q[,,1]),
					qs.net = as.vector(q[,,2]),
					qs.spear = as.vector(q[,,3]))
heatmap.df <- heatmap.df %>%
	mutate(sc = as.factor(sc)) %>%
	mutate(sc = relevel(sc, "[5-10)"))

heat.hook <- heatmap.df %>%
	dplyr::select(fg, sc, qs.hook) %>%
	ggplot() +
		geom_tile(aes(x = sc, y = fg, fill = qs.hook)) +
		scale_fill_viridis_c(name = expression(qs[line])) +
		labs(title = "Hook-and-line", x = "", y = "") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
			 axis.text.y = element_text(size = 8))

heat.net <- heatmap.df %>%
	dplyr::select(fg, sc, qs.net) %>%
	ggplot() +
		geom_tile(aes(x = sc, y = fg, fill = qs.net)) +
		scale_fill_viridis_c(name = expression(qs[net])) +
		labs(title = "Net", x = "", y = "") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
		theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
			 axis.text.y = element_text(size = 8))

heat.spear <- heatmap.df %>%
	dplyr::select(fg, sc, qs.spear) %>%
	ggplot() +
		geom_tile(aes(x = sc, y = fg, fill = qs.spear)) +
		scale_fill_viridis_c(name = expression(qs[spear])) +
		labs(title = "Spear", x = "", y = "") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
		theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
			 axis.text.y = element_text(size = 8))

ggarrange(heat.hook, heat.net, heat.spear, labels = c("A","B","C"), label.x = 0.06)


# --------------------------------------------------- FIG SX: functional group mean length ---------------------------------------------------

# out <- list(B.equil[1], cB.equil[2], Bfg.equil[3], cBfg.equil[4], BmeanW[5], BmeanL[6], BmeanLmat[7], 
#		    CmeanW[8],CmeanL[9], fg.BmeanW[10], fg.BmeanL[11], fg.BmeanLmat[12], fg.CmeanW[13], fg.CmeanL[14])

# Create dataframe with functional group biomass and catch for all scenarios
fgmeanL.df <- data.frame(scenario = rep(scenarios, each = (length(fg) * length(effort))),
					fg = rep(fg, each = length(effort), times = length(scenarios)),
					effort = rep(effort, times = (length(fg) * length(scenarios))),
					BmeanL = c(as.vector(t(scen.1.out[[11]][1,,])), as.vector(t(scen.2.out[[11]][1,,])), as.vector(t(scen.3.out[[11]][1,,])),
							 as.vector(t(scen.4.out[[11]][1,,])), as.vector(t(scen.5.out[[11]][1,,])), as.vector(t(scen.6.out[[11]][1,,])),
							 as.vector(t(scen.7.out[[11]][1,,]))),
					CmeanL = c(as.vector(t(scen.1.out[[14]][1,,])), as.vector(t(scen.2.out[[14]][1,,])), as.vector(t(scen.3.out[[14]][1,,])),
							 as.vector(t(scen.4.out[[14]][1,,])), as.vector(t(scen.5.out[[14]][1,,])), as.vector(t(scen.6.out[[14]][1,,])),
							 as.vector(t(scen.7.out[[14]][1,,]))))
					
# Iterate through functional groups to construct plots
fgmeanL.B <- vector(mode="list", length=length(fg)) # empty list to store ggplots
fgmeanL.C <- vector(mode="list", length=length(fg)) # empty list to store ggplots
for(i in 1:length(fg)){
   
   fg.i <- fg[i] # Save name of functional group
   
   # Create dataframe for fg.i
   fg.df <- fgmeanL.df %>% 
   	filter(fg == fg.i) %>%                       		 # get data for fg.i
   	mutate(effort.rel = effort/max(effort)) %>%  		 # relative effort
	mutate(BmeanL.rel = BmeanL/max(BmeanL, na.rm=TRUE)) %>% # relative biomass
   	mutate(CmeanL.rel = CmeanL/max(CmeanL, na.rm=TRUE))     # relative catch
   tmp.df <- plot1Half.df %>% 
   	left_join(., fg.df, by = c("effort.rel","scenario")) # get biomass and catch of fg.1 when total biomass = 0.5B0
   
   # create biomass plot for fg.i and all scenarios
   plot.fgmeanL.B <- ggplot() +
   	geom_line(data = fg.df, aes(x = effort.rel, y = BmeanL.rel, color = scenario), lwd = 0.5) +
   	geom_point(data = tmp.df, aes(x = effort.rel, y = BmeanL.rel, color = scenario), size = 1.5, alpha = 0.75) +
   	scale_y_continuous(expand = c(0.02,0.02)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
   	labs(x = "", y = "", title = fg.i) +
   	theme_classic() +
   	theme(plot.title = element_text(hjust = 0.5, size = 9),
   		 legend.title = element_blank(),
   		 legend.position = "none",
   		 plot.margin = margin(0.3,0.3,0.3,0.3,"cm"))
   	 
   fgmeanL.B[[i]] <- plot.fgmeanL.B
   
   # create catch plot for fg.i and all scenarios
   plot.fgmeanL.C <- ggplot() +
   	geom_line(data = fg.df, aes(x = effort.rel, y = CmeanL.rel, color = scenario), lwd = 0.5) +
   	geom_point(data = tmp.df, aes(x = effort.rel, y = CmeanL.rel, color = scenario), size = 1.5, alpha = 0.75) +
   	scale_y_continuous(expand = c(0.02,0.02)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
   	labs(x = "", y = "", title = fg.i) +
   	theme_classic() +
   	theme(plot.title = element_text(hjust = 0.5, size = 9),
   		 legend.title = element_blank(),
   		 legend.position = "none",
   		 plot.margin = margin(0.3,0.3,0.3,0.3,"cm"))
   
   fgmeanL.C[[i]] <- plot.fgmeanL.C
}

si.fgBmeanL <- ggarrange(fgmeanL.B[[1]],fgmeanL.B[[2]],fgmeanL.B[[3]],fgmeanL.B[[4]],fgmeanL.B[[5]],
	 			     fgmeanL.B[[6]],fgmeanL.B[[7]],fgmeanL.B[[8]],fgmeanL.B[[9]], nrow=3, ncol=3,
					labels=c("A","B","C","D","E","F","G","H","I"), common.legend = TRUE, legend = "bottom",
					label.x = 0.1)
annotate_figure(si.fgBmeanL,
			 bottom = text_grob("Effort", vjust=-6),
			 left = text_grob("Relative mean length (stock)", rot = 90, vjust = 1.5, hjust = 0.15))

si.fgCmeanL <- ggarrange(fgmeanL.C[[1]],fgmeanL.C[[2]],fgmeanL.C[[3]],fgmeanL.C[[4]],fgmeanL.C[[5]],
					fgmeanL.C[[6]],fgmeanL.C[[7]],fgmeanL.C[[8]],fgmeanL.C[[9]], nrow=3, ncol=3,
		     		labels=c("A","B","C","D","E","F","G","H","I"), common.legend = TRUE, legend = "bottom",
					label.x = 0.1)
annotate_figure(si.fgCmeanL,
			 bottom = text_grob("Effort", vjust=-6),
			 left = text_grob("Relative mean length (catch)", rot = 90, vjust = 1.5, hjust = 0.15))


# --------------------------------------------------- FIG SX: 0.5B0 sensitivity analysis ---------------------------------------------------

# Total Biomass 
totalB.df <- data.frame(Parameters = c("alpha", "beta", "mu", "sigma", "gamma"),
				    base = rep(base.sa$B, 5),
				    hi = c(alphaHi.sa$B, betaHi.sa$B, muHi.sa$B, sigmaHi.sa$B, GeHi.sa$B),
				    lo = c(alphaLo.sa$B, betaLo.sa$B, muLo.sa$B, sigmaLo.sa$B, GeLo.sa$B))
totalB.df$hi_change <- ((totalB.df$hi - totalB.df$base) / totalB.df$base) * 100
totalB.df$lo_change <- ((totalB.df$lo - totalB.df$base) / totalB.df$base) * 100

sa1 <- ggplot(data = totalB.df) +
  geom_errorbar(aes(x = Parameters, ymin = lo_change, ymax = hi_change, width = 0.25), lwd = 1) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray", alpha = 0.4) +
  labs(x = "", y = "", title = "Total biomass") +
  scale_x_discrete(labels = c(expression(alpha),expression(beta),expression(gamma),expression(mu),expression(sigma))) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

# Total Catch 
totalC.df <- data.frame(Parameters = c("alpha", "beta", "mu", "sigma", "gamma"),
				    base = rep(base.sa$cB, 5),
				    hi = c(alphaHi.sa$cB, betaHi.sa$cB, muHi.sa$cB, sigmaHi.sa$cB, GeHi.sa$cB),
				    lo = c(alphaLo.sa$cB, betaLo.sa$cB, muLo.sa$cB, sigmaLo.sa$cB, GeLo.sa$cB))
totalC.df$hi_change <- ((totalC.df$hi - totalC.df$base) / totalC.df$base) * 100
totalC.df$lo_change <- ((totalC.df$lo - totalC.df$base) / totalC.df$base) * 100

sa2 <- ggplot(data = totalC.df) +
  geom_errorbar(aes(x = Parameters, ymin = lo_change, ymax = hi_change, width = 0.25), lwd = 1) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray", alpha = 0.4) +
  labs(x = "", y = "", title = "Total catch") +
  scale_x_discrete(labels = c(expression(alpha),expression(beta),expression(gamma),expression(mu),expression(sigma))) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

# Mean stock length
BmeanL.df <- data.frame(Parameters = c("alpha", "beta", "mu", "sigma", "gamma"),
				    base = rep(base.sa$BmeanL, 5),
				    hi = c(alphaHi.sa$BmeanL, betaHi.sa$BmeanL, muHi.sa$BmeanL, sigmaHi.sa$BmeanL, GeHi.sa$BmeanL),
				    lo = c(alphaLo.sa$BmeanL, betaLo.sa$BmeanL, muLo.sa$BmeanL, sigmaLo.sa$BmeanL, GeLo.sa$BmeanL))
BmeanL.df$hi_change <- ((BmeanL.df$hi - BmeanL.df$base) / BmeanL.df$base) * 100
BmeanL.df$lo_change <- ((BmeanL.df$lo - BmeanL.df$base) / BmeanL.df$base) * 100

sa3 <- ggplot(data = BmeanL.df) +
  geom_errorbar(aes(x = Parameters, ymin = lo_change, ymax = hi_change, width = 0.25), lwd = 1) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray", alpha = 0.4) +
  labs(x = "", y = "", title = "Mean length (stock)") +
  scale_x_discrete(labels = c(expression(alpha),expression(beta),expression(gamma),expression(mu),expression(sigma))) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

# Mean catch length
CmeanL.df <- data.frame(Parameters = c("alpha", "beta", "mu", "sigma", "gamma"),
				    base = rep(base.sa$CmeanL, 5),
				    hi = c(alphaHi.sa$CmeanL, betaHi.sa$CmeanL, muHi.sa$CmeanL, sigmaHi.sa$CmeanL, GeHi.sa$CmeanL),
				    lo = c(alphaLo.sa$CmeanL, betaLo.sa$CmeanL, muLo.sa$CmeanL, sigmaLo.sa$CmeanL, GeLo.sa$CmeanL))
CmeanL.df$hi_change <- ((CmeanL.df$hi - CmeanL.df$base) / CmeanL.df$base) * 100
CmeanL.df$lo_change <- ((CmeanL.df$lo - CmeanL.df$base) / CmeanL.df$base) * 100

sa4 <- ggplot(data = CmeanL.df) +
  geom_errorbar(aes(x = Parameters, ymin = lo_change, ymax = hi_change, width = 0.25), lwd = 1) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray", alpha = 0.4) +
  labs(x = "", y = "", title = "Mean length (catch)") +
  scale_x_discrete(labels = c(expression(alpha),expression(beta),expression(gamma),expression(mu),expression(sigma))) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

sa.plot <- ggarrange(sa1, sa2, sa3, sa4, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"), label.x = 0.05)
annotate_figure(sa.plot,
			 bottom = text_grob("Parameters", vjust=0),
			 left = text_grob("Percent change", rot = 90, vjust = 1.5, hjust = 0.15))


##### TEST RECRUITMENT CHANGES

B <- seq(1,400)
R <- B / (alpha[1] + beta[1] * B)
plot(B, R, type = "l")
R1 <- B / (alpha[1] + (beta[1]+beta[1]*0.1) * B)
R2 <- B / (alpha[1] + (beta[1]-beta[1]*0.1) * B)
lines(B, R1, lty = "dashed", col = "green")
lines(B, R2, lty = "dotted", col = "red")
