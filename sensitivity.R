### sensitivity
###
### Author: Paul Carvalho
###
### Description: Run sensitivity analyses

sensitivity <- function(effort, total.effort, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q, sens, plots){
     # Model scenario 1
     gear.mgmt.1 <- c(1,1,1) # all gears used
     model.1 <- run_model(effort, gear.mgmt.1, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q)
     N.ijte.1 <- model.1[[1]]
     B.ijte.1 <- model.1[[2]]
     c.ijte.1 <- model.1[[3]]
     
     # Model scenario 2
     gear.mgmt.2 <- c(0,1,1) # no line
     model.2 <- run_model(effort, gear.mgmt.2, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q)
     N.ijte.2 <- model.2[[1]]
     B.ijte.2 <- model.2[[2]]
     c.ijte.2 <- model.2[[3]]
     
     # Model scenario 3
     gear.mgmt.3 <- c(1,0,1) # no net
     model.3 <- run_model(effort, gear.mgmt.3, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q)
     N.ijte.3 <- model.3[[1]]
     B.ijte.3 <- model.3[[2]]
     c.ijte.3 <- model.3[[3]]
     
     # Model scenario 4
     gear.mgmt.4 <- c(1,1,0) # no spear
     model.4 <- run_model(effort, gear.mgmt.4, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q)
     N.ijte.4 <- model.4[[1]]
     B.ijte.4 <- model.4[[2]]
     c.ijte.4 <- model.4[[3]]
     
     # Model scenario 5
     gear.mgmt.5 <- c(1,0,0) # only line
     model.5 <- run_model(effort, gear.mgmt.5, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q)
     N.ijte.5 <- model.5[[1]]
     B.ijte.5 <- model.5[[2]]
     c.ijte.5 <- model.5[[3]]
     
     # Model scenario 6
     gear.mgmt.6 <- c(0,1,0) # only net
     model.6 <- run_model(effort, gear.mgmt.6, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q)
     N.ijte.6 <- model.6[[1]]
     B.ijte.6 <- model.6[[2]]
     c.ijte.6 <- model.6[[3]]
     
     # Model scenario 7
     gear.mgmt.7 <- c(0,0,1) # only spear
     model.7 <- run_model(effort, gear.mgmt.7, nsc, nspecies, t, r, M, phi, L.lower, L.upper, W.a, W.b, q)
     N.ijte.7 <- model.7[[1]]
     B.ijte.7 <- model.7[[2]]
     c.ijte.7 <- model.7[[3]]
        
     # Summarize data
     summarize_data <- function(X, mgmt, Borc){
          tmp <- X[,,t,]
          tmp <- colSums(colSums(tmp, na.rm=TRUE))
          if(Borc == "B"){
               out <- data.frame(mgmt_sc = mgmt, B = mean(tmp, na.rm=TRUE))
          } else {
               out <- data.frame(mgmt_sc = mgmt, c = mean(tmp, na.rm=TRUE))
          }
          return(out)
     }
     
     B.df <- summarize_data(B.ijte.1, mgmt = "All gears", "B")
     B.df <- rbind(B.df,
                         summarize_data(B.ijte.2, mgmt = "Hook-and-line ban", "B"))
     B.df <- rbind(B.df,
                         summarize_data(B.ijte.3, mgmt = "Net ban", "B"))
     B.df <- rbind(B.df,
                         summarize_data(B.ijte.4, mgmt = "Spear ban", "B"))
     B.df <- rbind(B.df,
                         summarize_data(B.ijte.5, mgmt = "Hook-and-line", "B"))
     B.df <- rbind(B.df,
                         summarize_data(B.ijte.6, mgmt = "Net", "B"))
     B.df <- rbind(B.df,
                         summarize_data(B.ijte.7, mgmt = "Spear", "B"))
     
     c.df <- summarize_data(c.ijte.1, mgmt = "All gears", "c")
     c.df <- rbind(c.df,
                         summarize_data(c.ijte.2, mgmt = "Hook-and-line ban", "c"))
     c.df <- rbind(c.df,
                         summarize_data(c.ijte.3, mgmt = "Net ban", "c"))
     c.df <- rbind(c.df,
                         summarize_data(c.ijte.4, mgmt = "Spear ban", "c"))
     c.df <- rbind(c.df,
                         summarize_data(c.ijte.5, mgmt = "Hook-and-line", "c"))
     c.df <- rbind(c.df,
                         summarize_data(c.ijte.6, mgmt = "Net", "c"))
     c.df <- rbind(c.df,
                         summarize_data(c.ijte.7, mgmt = "Spear", "c"))
     
     names(B.df)[2] <- "V"
     B.df$Val <- "B"
     names(c.df)[2] <- "V"
     c.df$Val <- "c"
     
     out.df <- rbind(B.df,c.df)
     out.df$sensitivity = sens
     names(out.df)[1] <- "mgmt"
     
     # PLOTS ####
     if(plots == TRUE){
        output.B.1 <- calc_fgsize_output(B.ijte.1[,,t,], effort, total.effort, nsc, nspecies)
                B.rel.1 <- output.B.1[[1]]
                B.tot.1 <- output.B.1[[2]]
        output.B.2 <- calc_fgsize_output(B.ijte.2[,,t,], effort, total.effort, nsc, nspecies)
                B.rel.2 <- output.B.2[[1]]
                B.tot.2 <- output.B.2[[2]]                
        output.B.3 <- calc_fgsize_output(B.ijte.3[,,t,], effort, total.effort, nsc, nspecies)
                B.rel.3 <- output.B.3[[1]]
                B.tot.3 <- output.B.3[[2]]
        output.B.4 <- calc_fgsize_output(B.ijte.4[,,t,], effort, total.effort, nsc, nspecies)
                B.rel.4 <- output.B.4[[1]]
                B.tot.4 <- output.B.4[[2]]
        output.B.5 <- calc_fgsize_output(B.ijte.5[,,t,], effort, total.effort, nsc, nspecies)
                B.rel.5 <- output.B.5[[1]]
                B.tot.5 <- output.B.5[[2]]
        output.B.6 <- calc_fgsize_output(B.ijte.6[,,t,], effort, total.effort, nsc, nspecies)
                B.rel.6 <- output.B.6[[1]]
                B.tot.6 <- output.B.6[[2]]
        output.B.7 <- calc_fgsize_output(B.ijte.7[,,t,], effort, total.effort, nsc, nspecies)
                B.rel.7 <- output.B.7[[1]]
                B.tot.7 <- output.B.7[[2]]
        output.c.1 <- calc_fgsize_output(c.ijte.1[,,t,], effort, total.effort, nsc, nspecies)
                c.rel.1 <- output.c.1[[1]]
                c.tot.1 <- output.c.1[[2]]
        output.c.2 <- calc_fgsize_output(c.ijte.2[,,t,], effort, total.effort, nsc, nspecies)
                c.rel.2 <- output.c.2[[1]]
                c.tot.2 <- output.c.2[[2]]
        output.c.3 <- calc_fgsize_output(c.ijte.3[,,t,], effort, total.effort, nsc, nspecies)
                c.rel.3 <- output.c.3[[1]]
                c.tot.3 <- output.c.3[[2]]
        output.c.4 <- calc_fgsize_output(c.ijte.4[,,t,], effort, total.effort, nsc, nspecies)
                c.rel.4 <- output.c.4[[1]]
                c.tot.4 <- output.c.4[[2]]
        output.c.5 <- calc_fgsize_output(c.ijte.5[,,t,], effort, total.effort, nsc, nspecies)
                c.rel.5 <- output.c.5[[1]]
                c.tot.5 <- output.c.5[[2]]
        output.c.6 <- calc_fgsize_output(c.ijte.6[,,t,], effort, total.effort, nsc, nspecies)
                c.rel.6 <- output.c.6[[1]]
                c.tot.6 <- output.c.6[[2]]
        output.c.7 <- calc_fgsize_output(c.ijte.7[,,t,], effort, total.effort, nsc, nspecies)
                c.rel.7 <- output.c.7[[1]]
                c.tot.7 <- output.c.7[[2]]
                
        B.tot.1$mgmt <- "All gears"
        B.tot.2$mgmt <- "Hook-and-line ban"
        B.tot.3$mgmt <- "Net ban"
        B.tot.4$mgmt <- "Spear ban" 
        B.tot.5$mgmt <- "Hook-and-line" 
        B.tot.6$mgmt <- "Net" 
        B.tot.7$mgmt <- "Spear" 
        B.tot.all <- rbind(B.tot.1, B.tot.2, B.tot.3, B.tot.4, B.tot.5, B.tot.6, B.tot.7)
        B.tot.all$V <- B.tot.all$V / max(B.tot.all$V)
        all.plot.1 <- ggplot() +
                geom_line(data = B.tot.all, aes(x = E, y = V, color = mgmt), lwd = 1) +
                geom_hline(yintercept = 0.5, lty = "dashed", color = "gray10", lwd = 1, alpha = 0.5) +
                scale_color_manual(values = cbp2) +
                # geom_hline(yintercept = 0.25, lty = "dashed", color = "gray10", lwd = 1, alpha = 0.5) +
                # scale_x_continuous(limits = c(0,1.05), expand = c(0,0)) +
                # scale_y_continuous(limits = c(0,1.05), expand = c(0,0)) +
                theme_classic()
        
        # find effort that results in 0.5Bo
        B.tot.all$sqdiff <- ((B.tot.all$V-0.5)^2)
        B.half <- B.tot.all %>%
                group_by(mgmt) %>%
                summarise(min = min(sqdiff)) %>%
                left_join(., B.tot.all, by = c("min" = "sqdiff")) %>%
                dplyr::select(mgmt = mgmt.x, E = E, V = V)
         
        c.tot.1$mgmt <- "All gears"
        c.tot.2$mgmt <- "Hook-and-line ban"
        c.tot.3$mgmt <- "Net ban"
        c.tot.4$mgmt <- "Spear ban" 
        c.tot.5$mgmt <- "Hook-and-line" 
        c.tot.6$mgmt <- "Net" 
        c.tot.7$mgmt <- "Spear" 
        c.tot.all <- rbind(c.tot.1, c.tot.2, c.tot.3, c.tot.4, c.tot.5, c.tot.6, c.tot.7)
        c.tot.all$V <- c.tot.all$V / max(c.tot.all$V)

        B.half <- B.half %>%
                left_join(., c.tot.all, by = c("E","mgmt"))

        all.plot.2 <- ggplot() +
                geom_line(data = c.tot.all, aes(x = E, y = V, color = mgmt), lwd = 1.2, alpha = 0.5) +
                geom_point(data = B.half, aes(x = E, y = V.y, color = mgmt), size = 3, alpha = 0.5) +
                geom_linerange(data = B.half, aes(x = E, ymax = V.y, ymin = 0, color = mgmt), lty = "dashed", lwd = 1) +
                geom_linerangeh(data = B.half, aes(y = V.y, xmin = 0, xmax = E, color = mgmt), lty = "dashed", lwd = 1) +
                scale_color_manual(values = cbp2) +
                scale_y_continuous(expand = c(0,0), limits = c(0,1.05)) +
                scale_x_continuous(expand = c(0,0), limits = c(0,1.05)) +
                theme_classic() +
                theme(legend.position = "none")
        
        # PLOT BIOMASS AND CATCH FOR EACH FUNCTIONAL GROUP
        
        # Browsers
        brow.B.all <- fg_out_all("Browser", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        brow.c.all <- fg_out_all("Browser", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot1 <- fg_plot_all(brow.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot1 <- fg_plot_all(brow.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = TRUE)
        
        # Detritivores
        detr.B.all <- fg_out_all("Detritivore", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        detr.c.all <- fg_out_all("Detritivore", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot2 <- fg_plot_all(detr.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot2 <- fg_plot_all(detr.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        # Excavators/scrapers
        exsc.B.all <- fg_out_all("Excavator/scraper", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        exsc.c.all <- fg_out_all("Excavator/scraper", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot3 <- fg_plot_all(exsc.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot3 <- fg_plot_all(exsc.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        # Grazers
        graz.B.all <- fg_out_all("Grazer", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        graz.c.all <- fg_out_all("Grazer", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot4 <- fg_plot_all(graz.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot4 <- fg_plot_all(graz.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        # Macro-invertivore
        macro.B.all <- fg_out_all("Macro-invertivore", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        macro.c.all <- fg_out_all("Macro-invertivore", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot5 <- fg_plot_all(macro.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot5 <- fg_plot_all(macro.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        # Micro-invertivore
        micro.B.all <- fg_out_all("Micro-invertivore", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        micro.c.all <- fg_out_all("Micro-invertivore", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot6 <- fg_plot_all(micro.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot6 <- fg_plot_all(micro.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        # Pisci-invertivore
        pive.B.all <- fg_out_all("Pisci-invertivore", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        pive.c.all <- fg_out_all("Pisci-invertivore", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot7 <- fg_plot_all(pive.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot7 <- fg_plot_all(pive.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        # Piscivore
        pisc.B.all <- fg_out_all("Piscivore", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        pisc.c.all <- fg_out_all("Piscivore", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot8 <- fg_plot_all(pisc.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot8 <- fg_plot_all(pisc.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        # Planktivore
        plan.B.all <- fg_out_all("Planktivore", B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7) 
        plan.c.all <- fg_out_all("Planktivore", c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        B.plot9 <- fg_plot_all(plan.B.all, B.half, Borc = "B", x.lab = NULL, y.lab = NULL)
        c.plot9 <- fg_plot_all(plan.c.all, B.half, Borc = "c", x.lab = NULL, y.lab = NULL, legend = FALSE)
        
        
        ### SIZE DISTRIBUTIONS 
        B.sizedist.all <- sizedist_out_all(B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7)
        c.sizedist.all <- sizedist_out_all(c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7)
        sizedist.all <- B.sizedist.all %>%
                left_join(., c.sizedist.all, by = c("E", "mgmt")) %>%
                dplyr::select(E, mgmt, B = V.x, B.small = small.x, B.propsmall = prop.small.x,
                              c = V.y, c.small = small.y, c.propsmall = prop.small.y)

        # Plot
        b.propsmall.plot <- ggplot() +
                geom_line(data = sizedist.all, aes(x = E, y = B.propsmall, color = mgmt), lwd = 1) +
                labs(x = NULL, y = "Propotion of small fishes\n(biomass)") +
                scale_color_manual(values = cbp2) +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) +
                scale_x_continuous(limits = c(0, 1.05), expand = c(0,0)) +
                theme(legend.title = element_blank(),
                      legend.position = c(0.15, 0.75))
        c.propsmall.plot <- ggplot() +
                geom_line(data = sizedist.all, aes(x = E, y = c.propsmall, color = mgmt), lwd = 1) +
                labs(x = "Effort", y = "Propotion of small fishes\n(catch)") +
                scale_color_manual(values = cbp2) +
                theme_classic() +
                scale_y_continuous(expand = c(0,0)) +
                scale_x_continuous(limits = c(0, 1.05), expand = c(0,0)) +
                theme(legend.position = "none")
                
        ### SIZE DISTRIBUTION BY FUNCTIONAL GROUP
        
        B.sizedist.fg <- sizedist_out_all(B.rel.1, B.rel.2, B.rel.3, B.rel.4, B.rel.5, B.rel.6, B.rel.7, fg = TRUE)
        c.sizedist.fg <- sizedist_out_all(c.rel.1, c.rel.2, c.rel.3, c.rel.4, c.rel.5, c.rel.6, c.rel.7, fg = TRUE)
        sizedist.fg <- B.sizedist.fg %>%
                left_join(., c.sizedist.fg, by = c("mgmt","E","fg")) %>%
                dplyr::select(fg, E, mgmt, B = V.x, B.small = sizedist.x, c = V.y, c.small = sizedist.y)

        # Plot
        brow.size.plot.B <- fg_size_all("Browser", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = TRUE)
        detr.size.plot.B <- fg_size_all("Detritivore", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        exsc.size.plot.B <- fg_size_all("Excavator/scraper", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        graz.size.plot.B <- fg_size_all("Grazer", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        macro.size.plot.B <- fg_size_all("Macro-invertivore", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        micro.size.plot.B <- fg_size_all("Micro-invertivore", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        pive.size.plot.B <- fg_size_all("Pisci-invertivore", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        pisc.size.plot.B <- fg_size_all("Piscivore", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        plan.size.plot.B <- fg_size_all("Planktivore", sizedist.fg, "B", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        
        
        brow.size.plot.c <- fg_size_all("Browser", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = TRUE)
        detr.size.plot.c <- fg_size_all("Detritivore", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        exsc.size.plot.c <- fg_size_all("Excavator/scraper", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        graz.size.plot.c <- fg_size_all("Grazer", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        macro.size.plot.c <- fg_size_all("Macro-invertivore", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        micro.size.plot.c <- fg_size_all("Micro-invertivore", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        pive.size.plot.c <- fg_size_all("Pisci-invertivore", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        pisc.size.plot.c <- fg_size_all("Piscivore", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        plan.size.plot.c <- fg_size_all("Planktivore", sizedist.fg, "c", x.lab = NULL, 
                                      y.lab = "\n", legend = FALSE)
        
        
        out.df <- list(out.df, all.plot.1, all.plot.2,
                        B.plot1, B.plot2, B.plot3, B.plot4, B.plot5, B.plot6, B.plot7, B.plot8, B.plot9,
                        c.plot1, c.plot2, c.plot3, c.plot4, c.plot5, c.plot6, c.plot7, c.plot8, c.plot9,
                        b.propsmall.plot, c.propsmall.plot,
                        brow.size.plot.B, detr.size.plot.B, exsc.size.plot.B, graz.size.plot.B, macro.size.plot.B, 
                        micro.size.plot.B, pive.size.plot.B, pisc.size.plot.B, plan.size.plot.B,
                        brow.size.plot.c, detr.size.plot.c, exsc.size.plot.c, graz.size.plot.c, macro.size.plot.c, 
                        micro.size.plot.c, pive.size.plot.c, pisc.size.plot.c, plan.size.plot.c,
                        B.tot.all, B.half, c.tot.all)
        }
        return(out.df)
}


