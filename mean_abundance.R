### mean_abundance
###
### Author: Paul Carvalho
###
### Description: This script is called by "main_model.R" and calculates the mean
###              abundance per hectare for each functional group.

# Libraries ---------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)

# Functions ---------------------------------------------------------------
standard_error <- function(x){
     return(sd(x) / sqrt(length(x)))
}

# Load data ---------------------------------------------------------------
uvc.fish <- read.csv("indonesia_uvc_fish.csv")
fish.key <- read.csv("new_species_key.csv")

# remove unique values from fish.key
fish.key <- distinct(fish.key)

# assign temp functional groups
uvc.fish.1 <- uvc.fish %>% 
     dplyr::filter(!(site_name == "Sombano")) %>% # remove Sombano due to incomplete survey
     dplyr::select(site_name, region, transect, genus_species, genus, species, family, size_cm, size_5cm_bin, abundance, 
            old_fg = functional_group, a, b, observer) %>%
     mutate(genus_species = tolower(genus_species)) %>%
     left_join(., fish.key, by = "genus_species")

# refine functional groups
uvc.fish.2 <- uvc.fish.1 %>% 
        filter(fg != "Corallivore") %>% 
        filter(fg != "Spongivore") %>%
        mutate(fg = as.character(fg))

# Caclulate mean and standard error of abundance per hectare
fg.abundance <- uvc.fish.2 %>%
     filter(region == "wakatobi") %>%
     group_by(site_name, transect, fg) %>%
     summarise(abundance = sum(abundance)) %>%
     group_by(fg) %>%
     summarise(mean_ab = mean(abundance),
               sd_ab = sd(abundance),
               se_ab = standard_error(abundance)) %>%
     mutate(mean_ab = mean_ab * 40) %>% # convert from 250m2 to hectare
     mutate(sd_ab = sd_ab * 40) %>%
     mutate(se_ab = se_ab * 40)

# ggplot() +
#      geom_bar(data = fg.abundance, aes(x = reorder(func_group, -mean_ab), y = mean_ab), stat = "identity") +
#      geom_errorbar(data = fg.abundance, aes(x = reorder(func_group, -mean_ab), ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.3) +
#      labs(x = "Functional group", y = "Abundance per hectare") +
#      theme_classic() +
#      scale_y_continuous(expand = c(0,0))

# Calculate mean and standard error of biomass per hectare
fg.biomass <- uvc.fish.2 %>%
     filter(region == "wakatobi") %>%
     mutate(biomass_g = (a * size_cm ^ b) * abundance) %>%
     mutate(biomass_kg = biomass_g / 1000) %>%
     group_by(site_name, transect, fg) %>%
     summarise(biomass = sum(biomass_kg)) %>%
     group_by(fg) %>%
     summarise(mean_bio = mean(biomass),
               sd_bio = sd(biomass),
               se_bio = standard_error(biomass)) %>%
     mutate(mean_bio = mean_bio * 40) %>%
     mutate(sd_bio = sd_bio * 40) %>%
     mutate(se_bio = se_bio * 40)

# ggplot() +
#      geom_bar(data = fg.biomass, aes(x = reorder(func_group, -mean_bio), y = mean_bio), stat = "identity") +
#      geom_errorbar(data = fg.biomass, aes(x = reorder(func_group, -mean_bio), ymin = mean_bio - se_bio, ymax = mean_bio + se_bio), width = 0.3) +
#      labs(x = "Functional group", y = "Biomass (kg/ha)") +
#      theme_classic() +
#      scale_y_continuous(expand = c(0,0))

uvc.fish <- uvc.fish.2
rm(uvc.fish.1, uvc.fish.2, fish.key)
