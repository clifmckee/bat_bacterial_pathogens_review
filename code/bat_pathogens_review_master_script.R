##########################
# Code for determining which bat bacterial pathogen host species distributions overlap
# R version 4.2.1 "Funny-Looking Kid"
# Works: 2022-10-07
##########################

###################################
### Load packages and functions ###
###################################

# Necessary packages
library(here)
library(sp)
library(rgeos)
library(maps)
library(mapdata)
library(maptools)
library(geosphere)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggthemes)
library(reshape2)
library(RColorBrewer)
library(rgdal)
library(dplyr)
library(tidyr)
library(stringr)
library(mapproj)
library(broom)
library(raster)
library(rmapshaper)
library(readxl)

# Define functions
"%ni%" <- Negate("%in%")

# Check working directory
here()

###################
### Import data ###
###################

# Read in the shape file
# Available from IUCN at https://www.iucnredlist.org/resources/spatial-data-download (newest update 2022-07-21)
# readShapeSpatial() is in the "sp" package
# Shape files for the distributions of every terrestrial mammal in the IUCN database
# mammterr <- readOGR(dsn="./data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
# save(mammterr, file = "./data/MAMMALS_TERRESTRIAL_ONLY/mammterr.RData")
load(file = "./data/MAMMALS_TERRESTRIAL_ONLY/mammterr.RData")

###################
### Run scripts ###
###################

# Make map for all bat species
source("./code/bat_pathogens_review_bat_species.R")

# Make map for sampled bat species
source("./code/bat_pathogens_review_main.R")

# Combine plots
plotAB <- plot_grid(plotA, plotB, nrow = 2, rel_heights = c(0.35, 0.65))
# Save to file
ggsave(
  filename = "./results/map_overlap.pdf",
  plot = plotAB,
  device = "pdf",
  width = 16.67,
  height = 20,
  units = "in"
)
ggsave(
  filename = "./results/map_overlap.png",
  plot = plotAB,
  device = "png",
  dpi = 300,
  width = 16.67,
  height = 20,
  units = "in"
)
ggsave(
  filename = "./results/map_overlap.tiff",
  plot = plotAB,
  device = "tiff",
  dpi = 300,
  width = 16.67,
  height = 20,
  units = "in"
)

###################
### End of code ###
###################
