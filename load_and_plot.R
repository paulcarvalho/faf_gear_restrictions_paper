# ------------------------------------------------------------------------------------------------------------------
## File 'load_and_plot.R'
##
## Publication: Carvalho, PG and Humphries AT. (2021). Gear restrictions create conservation and fisheries trade-offs for 
##              management. Fish and Fisheries. DOI:https://doi.org/10.1111/faf.12607 
##
## Description: Create plots presented in the paper. Coral reef fisheries model of nine functional groups to test biomass
##              and catch response to various gear-based management scenarios. The model losely represents the Wakatobi 
##              coral reef fishery in Indonesia, where certain model parameters were derived.

# Clean workspace
rm(list = ls())

# --------------------------------------------------- LIBRARIES ---------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pracma)
library(openxlsx)
library(DescTools)

# --------------------------------------------------- LOAD MODEL DATA ---------------------------------------------------
# The workspace file is too large for uploading to github, but the file will be provided upon request until until we can find a way to compress the file and make it available in the repository.
load("model_data.RData")

# --------------------------------------------------- FIG 4: total biomass and catch ---------------------------------------------------

# Colorblind friendly palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# Create dataframe for total biomass and catch of all scenarios
scenarios <- c("All gears", "Hook-and-line ban", "Net ban", "Spear ban", "Hook-and-line", "Net", "Spear")
plot1.df  <- data.frame(scenario = rep(scenarios, each = length(effort)),
				    effort   = rep(effort, length(scenarios)),
				    biomass  = c(scen.1.out[[1]], scen.2.out[[1]], scen.3.out[[1]], scen.4.out[[1]], scen.5.out[[1]], scen.6.out[[1]], scen.7.out[[1]]),
				    catch    = c(scen.1.out[[2]], scen.2.out[[2]], scen.3.out[[2]], scen.4.out[[2]], scen.5.out[[2]], scen.6.out[[2]], scen.7.out[[2]]))
bs.index <- c(6, 10, 15)
plot1.iqr <- data.frame(scenario = rep(scenarios, each = 3),
                       effort   = rep(effort.1.bs, length(scenarios)),
                       B.mu     = c(scen.1.out[[1]][bs.index], scen.2.out[[1]][bs.index], scen.3.out[[1]][bs.index], scen.4.out[[1]][bs.index], scen.5.out[[1]][bs.index], scen.6.out[[1]][bs.index], scen.7.out[[1]][bs.index]),
                       CB.mu    = c(scen.1.out[[2]][bs.index], scen.2.out[[2]][bs.index], scen.3.out[[2]][bs.index], scen.4.out[[2]][bs.index], scen.5.out[[2]][bs.index], scen.6.out[[2]][bs.index], scen.7.out[[2]][bs.index]),
                       B.iqr    = c(scen.1.bs[[1]]$B.iqr, scen.2.bs[[1]]$B.iqr, scen.3.bs[[1]]$B.iqr, scen.4.bs[[1]]$B.iqr, scen.5.bs[[1]]$B.iqr, scen.6.bs[[1]]$B.iqr, scen.7.bs[[1]]$B.iqr),
                       CB.iqr   = c(scen.1.bs[[1]]$CB.iqr, scen.2.bs[[1]]$CB.iqr, scen.3.bs[[1]]$CB.iqr, scen.4.bs[[1]]$CB.iqr, scen.5.bs[[1]]$CB.iqr, scen.6.bs[[1]]$CB.iqr, scen.7.bs[[1]]$CB.iqr))
                       
# Convert to relative values
plot1.df <- plot1.df %>%
	dplyr::mutate(effort.rel = effort/max(effort)) %>%    # relative effort
     dplyr::mutate(biomass.rel = biomass/max(biomass)) %>% # relative biomass
     dplyr::mutate(catch.rel = catch/max(catch)) %>%       # relative catch
     dplyr::mutate(sq.diff = (biomass.rel-0.5)^2)          # find squared difference between relative biomass and 0.5 to identify effort and catch at 0.5B0
plot1.iqr <- plot1.iqr %>%
     dplyr::mutate(effort.rel = effort/max(plot1.df$effort)) %>% # relative effort
     dplyr::mutate(B.mu = B.mu/max(plot1.df$biomass)) %>%
     dplyr::mutate(B.iqr = B.iqr/max(plot1.df$biomass)) %>%
     dplyr::mutate(CB.mu = CB.mu/max(plot1.df$catch)) %>%
     dplyr::mutate(CB.iqr = CB.iqr/max(plot1.df$catch))
     
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
     geom_point(aes(x = plot1.iqr$effort.rel, y = plot1.iqr$B.mu, color = plot1.iqr$scenario), size = 0.75) +     
     geom_errorbar(aes(x = plot1.iqr$effort.rel, ymin = plot1.iqr$B.mu - (plot1.iqr$B.iqr/2), ymax = plot1.iqr$B.mu + (plot1.iqr$B.iqr/2), color = plot1.iqr$scenario), lwd = 0.5, width = 0.01) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
	labs(x = "", y = expression(paste("Stock biomass (relative to ", B[paste(" ", italic("F"),"=0")], ")"))) +
	theme_classic() +
	theme(legend.title = element_blank(), legend.position = c(0.85,0.8), plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))

plot1.C <- ggplot() +
	geom_line(data = plot1.df, aes(x = effort.rel, y = catch.rel, color = scenario), lwd = 0.5) +
     geom_point(aes(x = plot1.iqr$effort.rel, y = plot1.iqr$CB.mu, color = plot1.iqr$scenario), size = 0.75) +
     geom_errorbar(aes(x = plot1.iqr$effort.rel, ymin = plot1.iqr$CB.mu - (plot1.iqr$CB.iqr/2), ymax = plot1.iqr$CB.mu + (plot1.iqr$CB.iqr/2), color = plot1.iqr$scenario), lwd = 0.5, width = 0.01) +
     scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
	labs(x = "Effort", y = expression(paste("Catch (relative to ", C[italic("max")], " across all scenarios)"))) +
	theme_classic() +
	theme(legend.position = "none", plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
	
ggarrange(plot1.B, plot1.C, nrow=2, ncol=1, labels=c("A","B"))

# --------------------------------------------------- ANOVA 1 - total biomass and catch ---------------------------------------------------

# Create dataframe
anova1.df <- NULL
for(i in 1:length(scenarios)){
     tmp.df <- get(paste('scen.', i, '.bs', sep = ''))[[3]] %>% dplyr::select('effort', 'Btot', 'CBtot') %>% dplyr::mutate(scenario = scenarios[i]) %>% dplyr::mutate(effort = as.factor(effort), scenario = as.factor(scenario)) 
     anova1.df <- rbind(anova1.df, tmp.df)
}
anova1.Bres <- aov(Btot ~ effort * scenario, data = anova1.df)
EtaSq(anova1.Bres, anova = TRUE)

anova1.CBres <- aov(CBtot ~ effort * scenario, data = anova1.df)
EtaSq(anova1.CBres, anova = TRUE)

# --------------------------------------------------- FIG 3-4: functional group biomass and catch ---------------------------------------------------

# Create dataframe with functional group biomass and catch for all scenarios
fg <- c("Browsers", "Detritivores", "Excavators/scrapers", "Grazers", "Macro-invertivores", "Micro-invertivores", "Pisci-invertivores", "Piscivores", "Planktivores") 
plot3.df <- data.frame(scenario = rep(scenarios, each = (length(fg) * length(effort))),
				   fg = rep(fg, each = length(effort), times = length(scenarios)),
				   effort = rep(effort, times = (length(fg) * length(scenarios))),
				   biomass = c(as.vector(t(scen.1.out[[3]])), as.vector(t(scen.2.out[[3]])), as.vector(t(scen.3.out[[3]])), as.vector(t(scen.4.out[[3]])), as.vector(t(scen.5.out[[3]])), as.vector(t(scen.6.out[[3]])), as.vector(t(scen.7.out[[3]]))),
				   catch = c(as.vector(t(scen.1.out[[4]])), as.vector(t(scen.2.out[[4]])), as.vector(t(scen.3.out[[4]])), as.vector(t(scen.4.out[[4]])), as.vector(t(scen.5.out[[4]])), as.vector(t(scen.6.out[[4]])), as.vector(t(scen.7.out[[4]]))))
plot3.iqr <- data.frame(scenario = rep(scenarios, each = (length(fg) * length(bs.index))),
                       fg = rep(fg, each = length(bs.index), times = length(scenarios)),
                       effort = rep(effort[bs.index], times = (length(fg) * length(scenarios))),
                       biomass = c(as.vector(t(scen.1.out[[3]][,bs.index])), as.vector(t(scen.2.out[[3]][,bs.index])), as.vector(t(scen.3.out[[3]][,bs.index])), as.vector(t(scen.4.out[[3]][,bs.index])), as.vector(t(scen.5.out[[3]][,bs.index])), as.vector(t(scen.6.out[[3]][,bs.index])), as.vector(t(scen.7.out[[3]][,bs.index]))),
                       catch = c(as.vector(t(scen.1.out[[4]][,bs.index])), as.vector(t(scen.2.out[[4]][,bs.index])), as.vector(t(scen.3.out[[4]][,bs.index])), as.vector(t(scen.4.out[[4]][,bs.index])), as.vector(t(scen.5.out[[4]][,bs.index])), as.vector(t(scen.6.out[[4]][,bs.index])), as.vector(t(scen.7.out[[4]][,bs.index]))),
                       B.iqr = c(scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$B.iqr, scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$B.iqr, scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$B.iqr, scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$B.iqr, scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$B.iqr, scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$B.iqr, scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$B.iqr),
                       CB.iqr = c(scen.1.bs[[2]][order(scen.1.bs[[2]]$fg),]$CB.iqr, scen.2.bs[[2]][order(scen.2.bs[[2]]$fg),]$CB.iqr, scen.3.bs[[2]][order(scen.3.bs[[2]]$fg),]$CB.iqr, scen.4.bs[[2]][order(scen.4.bs[[2]]$fg),]$CB.iqr, scen.5.bs[[2]][order(scen.5.bs[[2]]$fg),]$CB.iqr, scen.6.bs[[2]][order(scen.6.bs[[2]]$fg),]$CB.iqr, scen.7.bs[[2]][order(scen.7.bs[[2]]$fg),]$CB.iqr))
                       
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
   fg.iqr <- plot3.iqr %>%
     filter(fg == fg.i) %>%
     mutate(effort.rel = effort/max(fg.df$effort)) %>%
     mutate(biomass.rel = biomass/max(fg.df$biomass)) %>%
     mutate(catch.rel = catch/max(fg.df$catch)) %>%
     mutate(B.iqr.lo = (biomass - (B.iqr/2))/max(fg.df$biomass)) %>%
     mutate(B.iqr.hi = (biomass + (B.iqr/2))/max(fg.df$biomass)) %>%
     mutate(CB.iqr.lo = (catch - (CB.iqr/2))/max(fg.df$catch)) %>%
     mutate(CB.iqr.hi = (catch + (CB.iqr/2))/max(fg.df$catch))
     
   # force lower limits to zero
   fg.iqr$B.iqr.lo[which(fg.iqr$B.iqr.lo < 0)] <- 0
   fg.iqr$CB.iqr.lo[which(fg.iqr$CB.iqr.lo < 0)] <- 0
   
   tmp.df <- plot1Half.df %>%
   	left_join(., fg.df, by = c("effort.rel","scenario")) # get biomass and catch of fg.1 when total biomass = 0.5B0
   
   # create biomass plot for fg.i and all scenarios
   fg.plotB <- ggplot() +
   	geom_line(data = fg.df, aes(x = effort.rel, y = biomass.rel, color = scenario), lwd = 0.35) +
        # geom_point(data = fg.df, aes(x = effort.rel, y = biomass.rel, color = scenario)) +
     geom_point(data = fg.iqr, aes(x = effort.rel, y = biomass.rel, color = scenario), size = 0.25) +
     geom_errorbar(data = fg.iqr, aes(x = effort.rel, ymin = B.iqr.lo, ymax = B.iqr.hi, color = scenario), lwd = 0.35, width = 0.02, alpha = 0.5) +
   	scale_y_continuous(expand = c(0,0), limits = c(0, 1.02)) +
	scale_x_continuous(expand = c(0,0)) +
	scale_color_manual(values = cbp2) +
   	labs(x = "", y = "", title = fg.i) +
   	theme_classic() +
   	theme(plot.title = element_text(hjust = 0.5, size = 9), legend.title = element_blank(), legend.position = "none", plot.margin = margin(0.4,0.4,0.4,0.4,"cm"))
   	 
   plot3.B[[i]] <- fg.plotB
   
   # create catch plot for fg.i and all scenarios
   if(max(fg.iqr$CB.iqr.hi) > 1) iqr.lim <- max(fg.iqr$CB.iqr.hi) else iqr.lim <- 1
   
   fg.plotC <- ggplot() +
   	geom_line(data = fg.df, aes(x = effort.rel, y = catch.rel, color = scenario), lwd = 0.35) +
     geom_point(data = fg.iqr, aes(x = effort.rel, y = catch.rel, color = scenario), size = 0.25) +
     geom_errorbar(data = fg.iqr, aes(x = effort.rel, ymin = CB.iqr.lo, ymax = CB.iqr.hi, color = scenario), lwd = 0.35, width = 0.02) +
   	scale_y_continuous(expand = c(0,0), limits = c(0, iqr.lim)) +
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

# --------------------------------------------------- NEW PROPORTIONAL LOSS FIGURE ---------------------------------------------------
cbp3 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255")

tmp.line.plot <- plot3.df %>% 
                 dplyr::filter(scenario == 'Hook-and-line') %>%
                 dplyr::group_by(effort) %>%
                 dplyr::mutate(B.prop = biomass / sum(biomass)) %>%
                 dplyr::mutate(C.prop = catch / sum(catch)) %>%
                 dplyr::mutate(B.cumsum = cumsum(B.prop))
tmp.line.plot$fg <- factor(tmp.line.plot$fg, levels = c("Planktivores", "Piscivores", "Pisci-invertivores", "Micro-invertivores", "Macro-invertivores", "Grazers", "Excavators/scrapers", "Detritivores",  "Browsers") )
prop.line.plot <- ggplot() +
     geom_area(aes(x = tmp.line.plot$effort/max(tmp.line.plot$effort), y = tmp.line.plot$B.prop, fill = tmp.line.plot$fg)) +
     scale_fill_manual(values = cbp3) +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     labs(title = 'Hook-and-line', x = '', y = 'Proportion of total stock biomass') +
     theme_classic()

tmp.net.plot <- plot3.df %>% 
     dplyr::filter(scenario == 'Net') %>%
     dplyr::group_by(effort) %>%
     dplyr::mutate(B.prop = biomass / sum(biomass)) %>%
     dplyr::mutate(C.prop = catch / sum(catch)) %>%
     dplyr::mutate(B.cumsum = cumsum(B.prop))
tmp.net.plot$fg <- factor(tmp.net.plot$fg, levels = c("Planktivores", "Piscivores", "Pisci-invertivores", "Micro-invertivores", "Macro-invertivores", "Grazers", "Excavators/scrapers", "Detritivores",  "Browsers") )
prop.net.plot <- ggplot() +
     geom_area(aes(x = tmp.net.plot$effort/max(tmp.net.plot$effort), y = tmp.net.plot$B.prop, fill = tmp.net.plot$fg)) +
     scale_fill_manual(values = cbp3) +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     labs(title = 'Net', x = 'Effort', y = '') +
     theme_classic()

tmp.spear.plot <- plot3.df %>% 
     dplyr::filter(scenario == 'Spear') %>%
     dplyr::group_by(effort) %>%
     dplyr::mutate(B.prop = biomass / sum(biomass)) %>%
     dplyr::mutate(C.prop = catch / sum(catch)) %>%
     dplyr::mutate(B.cumsum = cumsum(B.prop))
tmp.spear.plot$fg <- factor(tmp.spear.plot$fg, levels = c("Planktivores", "Piscivores", "Pisci-invertivores", "Micro-invertivores", "Macro-invertivores", "Grazers", "Excavators/scrapers", "Detritivores",  "Browsers") )
prop.spear.plot <- ggplot() +
     geom_area(aes(x = tmp.spear.plot$effort/max(tmp.spear.plot$effort), y = tmp.spear.plot$B.prop, fill = tmp.spear.plot$fg)) +
     scale_fill_manual(values = cbp3) +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     labs(title = 'Spear', x = 'Effort', y = 'Proportion of total stock biomass') +
     theme_classic() +
     theme(legend.position = 'right', legend.title = element_blank())

leg <- get_legend(prop.spear.plot)
prop.line.plot <- prop.line.plot + theme(legend.position = 'none', plot.margin=unit(c(0.1, 0.3, 0, 0.2),"cm"))
prop.net.plot <- prop.net.plot + theme(legend.position = 'none', plot.margin=unit(c(0.1, 0.3, 0, 0.2),"cm"))
prop.spear.plot <- prop.spear.plot + theme(legend.position = 'none', plot.margin=unit(c(0.1, 0.3, 0, 0.2),"cm"))

plot.prop <- ggarrange(prop.line.plot, prop.net.plot, prop.spear.plot, leg, nrow=2, ncol=2,
                   labels=c("A","B","C"))

tmp.net.plot %>% dplyr::filter(fg == 'Planktivores')

cprop.line <- ggplot() +
     geom_area(aes(x = tmp.line.plot$effort/max(tmp.line.plot$effort), y = tmp.line.plot$C.prop, fill = tmp.line.plot$fg)) +
     scale_fill_manual(values = cbp3) +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     labs(title = 'Hook-and-line', x = '', y = 'Proportion of total catch') +
     theme_classic() +
     theme(legend.position = 'none')
cprop.net <- ggplot() +
     geom_area(aes(x = tmp.net.plot$effort/max(tmp.net.plot$effort), y = tmp.net.plot$C.prop, fill = tmp.net.plot$fg)) +
     scale_fill_manual(values = cbp3) +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     labs(title = 'Net', x = '', y = 'Proportion of total catch') +
     theme_classic() +
     theme(legend.position = 'none')
cprop.spear <- ggplot() +
     geom_area(aes(x = tmp.spear.plot$effort/max(tmp.spear.plot$effort), y = tmp.spear.plot$C.prop, fill = tmp.spear.plot$fg)) +
     scale_fill_manual(values = cbp3) +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_continuous(expand = c(0, 0)) +
     labs(title = 'Spear', x = '', y = 'Proportion of total catch') +
     theme_classic() + 
     theme(legend.position = 'none')     
plot.cprop <- ggarrange(cprop.line, cprop.net, cprop.spear, leg, nrow=2, ncol=2,
                       labels=c("A","B","C"))



# --------------------------------------------------- ANOVA 2 - functional group biomass and catch ---------------------------------------------------

# Create dataframe for each functional group
fg.eta.B  <- vector(mode = 'list', length = length(fg))
fg.eta.CB <- vector(mode = 'list', length = length(fg))
fg_eta <- createWorkbook() 

for(j in 1:length(fg)){
     tmp.aov.df <- NULL
     for(i in 1:length(scenarios)){
          tmp.df <- get(paste('scen.', i, '.bs', sep = ''))[[4]] %>% dplyr::filter(fg == fg[j]) %>% dplyr::select('effort', 'fg', 'B', 'CB') %>% dplyr::mutate(scenario = scenarios[i]) %>% dplyr::mutate(effort = as.factor(effort), scenario = as.factor(scenario))
          tmp.aov.df <- rbind(tmp.aov.df, tmp.df)
     }
     tmp.aov.Bres <- aov(B ~ effort + scenario, data = tmp.aov.df)
     tmp.eta.Bres <- EtaSq(tmp.aov.Bres, anova = TRUE)
     
     tmp.aov.CBres <- aov(CB ~ effort + scenario, data = tmp.aov.df)
     tmp.eta.CBres <- EtaSq(tmp.aov.CBres, anova = TRUE)
     
     fg.eta.B[[j]]  <- tmp.eta.Bres
     fg.eta.CB[[j]] <- tmp.eta.CBres
     
     names(fg.eta.B)[j]  <- fg[j]
     names(fg.eta.CB)[j] <- fg[j]
     
     # estimated marginal means
     # lm1.B <- lm(B ~ effort * scenario, data = tmp.aov.df)
     # emm1 <- emmeans::emmeans(lm1.B, specs = "effort", by = "scenario")
     # emm1.df <- as.data.frame(emm1)
     # plot(emm1)
     #
     
     addWorksheet(fg_eta, fg[j])
     writeData(fg_eta, fg[j], cbind(tmp.eta.Bres, tmp.eta.CBres))
}
saveWorkbook(fg_eta, 'fg_eta.xlsx')

fg_table <- createWorkbook() 
for(i in 1:length(fg)){
     fg.i <- fg[i] # Save name of functional group
     
     # Create dataframe for fg.i
     fg.df <- plot3.df %>% 
          filter(fg == fg.i) %>%                         # get data for fg.i
          mutate(effort.rel = effort/max(effort)) %>%    # relative effort
          mutate(biomass.rel = biomass/max(biomass)) %>% # relative biomass
          mutate(catch.rel = catch/max(catch))           # relative catch
     
     # Create dataframe for IQR
     fg.iqr <- plot3.iqr %>%
          filter(fg == fg.i) %>%
          mutate(effort.rel = effort/max(fg.df$effort)) %>%
          mutate(biomass.rel = biomass/max(fg.df$biomass)) %>%
          mutate(catch.rel = catch/max(fg.df$catch)) %>%
          mutate(B.iqr.lo = (biomass - (B.iqr/2))/max(fg.df$biomass)) %>%
          mutate(B.iqr.hi = (biomass + (B.iqr/2))/max(fg.df$biomass)) %>%
          mutate(CB.iqr.lo = (catch - (CB.iqr/2))/max(fg.df$catch)) %>%
          mutate(CB.iqr.hi = (catch + (CB.iqr/2))/max(fg.df$catch))
     
     # force lower limits to zero
     fg.iqr$B.iqr.lo[which(fg.iqr$B.iqr.lo < 0)] <- 0
     fg.iqr$CB.iqr.lo[which(fg.iqr$CB.iqr.lo < 0)] <- 0
     fg.iqr$B.iqr <- paste0(round(fg.iqr$B.iqr.lo, 3), '-', round(fg.iqr$B.iqr.hi, 3))
     fg.iqr$CB.iqr <- paste0(round(fg.iqr$CB.iqr.lo, 3), '-', round(fg.iqr$CB.iqr.hi, 3))
     
     # remove unwanted variables
     fg.iqr <- fg.iqr %>% dplyr::select(scenario, effort, biomass.rel, B.iqr, catch.rel, CB.iqr) %>% dplyr::mutate(biomass.rel = round(biomass.rel, 3), catch.rel = round(catch.rel, 3))
     fg.iqr <- reshape2::dcast(reshape2::melt(fg.iqr, id.vars = c('scenario', 'effort')), scenario ~ effort + variable)
     fg.iqr$fg <- fg[i]
     fg.iqr <- fg.iqr[, c(14, 1:13)]
     
     addWorksheet(fg_table, fg[i])
     writeData(fg_table, fg[i], fg.iqr)
}
saveWorkbook(fg_table, 'irq_table.xlsx')

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
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
			 axis.text.y = element_text(size = 5),
			 plot.title = element_text(size = 10),
			 legend.title = element_text(size = 8),
			 legend.key.width = unit(0.5, 'cm'),
			 legend.key.height = unit(0.5, 'cm'))

heat.net <- heatmap.df %>%
	dplyr::select(fg, sc, qs.net) %>%
	ggplot() +
		geom_tile(aes(x = sc, y = fg, fill = qs.net)) +
		scale_fill_viridis_c(name = expression(qs[net])) +
		labs(title = "Net", x = "Size class (cm)", y = "") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
			 axis.text.y = element_text(size = 5),
			 plot.title = element_text(size = 10),
			 legend.title = element_text(size = 8),
			 legend.key.width = unit(0.5, 'cm'),
			 legend.key.height = unit(0.5, 'cm'))

heat.spear <- heatmap.df %>%
	dplyr::select(fg, sc, qs.spear) %>%
	ggplot() +
		geom_tile(aes(x = sc, y = fg, fill = qs.spear)) +
		scale_fill_viridis_c(name = expression(qs[spear])) +
		labs(title = "Spear", x = "Size class (cm)", y = "") +
		scale_x_discrete(expand = c(0,0)) +
		scale_y_discrete(expand = c(0,0)) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
			 axis.text.y = element_text(size = 5),
			 plot.title = element_text(size = 10),
			 legend.title = element_text(size = 8),
			 legend.key.width = unit(0.5, 'cm'),
			 legend.key.height = unit(0.5, 'cm'))

ggarrange(heat.hook, heat.net, heat.spear, labels = c("A","B","C"), label.x = 0)


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
  labs(x = "", y = "Percent change", title = "Total biomass") +
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
  labs(x = "", y = "Percent change", title = "Total catch") +
  scale_x_discrete(labels = c(expression(alpha),expression(beta),expression(gamma),expression(mu),expression(sigma))) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

# # Mean stock length
# BmeanL.df <- data.frame(Parameters = c("alpha", "beta", "mu", "sigma", "gamma"),
# 				    base = rep(base.sa$BmeanL, 5),
# 				    hi = c(alphaHi.sa$BmeanL, betaHi.sa$BmeanL, muHi.sa$BmeanL, sigmaHi.sa$BmeanL, GeHi.sa$BmeanL),
# 				    lo = c(alphaLo.sa$BmeanL, betaLo.sa$BmeanL, muLo.sa$BmeanL, sigmaLo.sa$BmeanL, GeLo.sa$BmeanL))
# BmeanL.df$hi_change <- ((BmeanL.df$hi - BmeanL.df$base) / BmeanL.df$base) * 100
# BmeanL.df$lo_change <- ((BmeanL.df$lo - BmeanL.df$base) / BmeanL.df$base) * 100
# 
# sa3 <- ggplot(data = BmeanL.df) +
#   geom_errorbar(aes(x = Parameters, ymin = lo_change, ymax = hi_change, width = 0.25), lwd = 1) +
#   theme_classic() +
#   geom_hline(yintercept = 0, color = "gray", alpha = 0.4) +
#   labs(x = "", y = "", title = "Mean length (stock)") +
#   scale_x_discrete(labels = c(expression(alpha),expression(beta),expression(gamma),expression(mu),expression(sigma))) +
#   theme(axis.text = element_text(size = 10),
#         axis.title = element_text(size = 10))
# 
# # Mean catch length
# CmeanL.df <- data.frame(Parameters = c("alpha", "beta", "mu", "sigma", "gamma"),
# 				    base = rep(base.sa$CmeanL, 5),
# 				    hi = c(alphaHi.sa$CmeanL, betaHi.sa$CmeanL, muHi.sa$CmeanL, sigmaHi.sa$CmeanL, GeHi.sa$CmeanL),
# 				    lo = c(alphaLo.sa$CmeanL, betaLo.sa$CmeanL, muLo.sa$CmeanL, sigmaLo.sa$CmeanL, GeLo.sa$CmeanL))
# CmeanL.df$hi_change <- ((CmeanL.df$hi - CmeanL.df$base) / CmeanL.df$base) * 100
# CmeanL.df$lo_change <- ((CmeanL.df$lo - CmeanL.df$base) / CmeanL.df$base) * 100
# 
# sa4 <- ggplot(data = CmeanL.df) +
#   geom_errorbar(aes(x = Parameters, ymin = lo_change, ymax = hi_change, width = 0.25), lwd = 1) +
#   theme_classic() +
#   geom_hline(yintercept = 0, color = "gray", alpha = 0.4) +
#   labs(x = "", y = "", title = "Mean length (catch)") +
#   scale_x_discrete(labels = c(expression(alpha),expression(beta),expression(gamma),expression(mu),expression(sigma))) +
#   theme(axis.text = element_text(size = 10),
#         axis.title = element_text(size = 10))

sa.plot <- ggarrange(sa1, sa2, nrow = 2, ncol = 1, labels = c("A", "B"))
annotate_figure(sa.plot,
			 bottom = text_grob("Parameters", vjust=0))


##### TEST RECRUITMENT CHANGES

B <- seq(1,400)
R <- B / (alpha[1] + beta[1] * B)
plot(B, R, type = "l")
R1 <- B / (alpha[1] + (beta[1]+beta[1]*0.1) * B)
R2 <- B / (alpha[1] + (beta[1]-beta[1]*0.1) * B)
lines(B, R1, lty = "dashed", col = "green")
lines(B, R2, lty = "dotted", col = "red")
