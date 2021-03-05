### calc_fgsize_output
###
### Author: Paul Carvalho
###
### Description: Calculate the biomass or catch and size distribution for each functional group
###              and effort combination

calc_fgsize_output <- function(X, effort, total.effort, nsc, nspecies){
        # empty dataframes for saving data
        out.df <- data.frame(fg = as.character(), E = as.numeric(), V = as.numeric(), V.rel = as.numeric(), sizedist = as.numeric())
        out.df2 <- data.frame(E = as.numeric(), V = as.numeric())
        
        func.groups <- c("Browser","Detritivore","Excavator/scraper","Grazer","Macro-invertivore",
                         "Micro-invertivore","Pisci-invertivore","Piscivore","Planktivore")
        
        for(i in 1:length(effort)){
                effort.tmp <- rep((effort[i] / total.effort), nspecies)
                X.tmp <- X[,,i]
                
                # calculate proportions in each size class
                X.sum <- repmat(t(as.matrix(colSums(X.tmp))), nsc, 1)
                X.prop <- X.tmp / X.sum
                # get only proportion in the smallest two size classes, that is
                # the size class of recruits and smallest size class exposed to fishing
                X.prop <- colSums(X.prop[1:2,])
                
                # calculate total biomass
                V.tmp <- colSums(X.tmp)
                V2.tmp <- sum(V.tmp)
                
                # calculate relative biomass for each functional group
                V.rel.tmp <- V.tmp / sum(V.tmp)
                
                # save all data
                out.tmp <- data.frame(fg=func.groups, E = effort.tmp, V = V.tmp, V.rel = V.rel.tmp, sizedist = X.prop)
                out.tmp2 <- data.frame(E = effort.tmp[1], V = V2.tmp)
                # combine with other data
                out.df <- rbind(out.df, out.tmp)
                out.df2 <- rbind(out.df2, out.tmp2)
        }
        
        return(list(out.df, out.df2))
}