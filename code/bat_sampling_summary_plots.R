##########################
# Code for plotting bacterial pathogen biogeography and sampling patterns among bats
# R version 4.2.1 "Funny-Looking Kid"
# Works: 2022-11-07
##########################

# # Installing ggtree package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtree")

# Necessary packages
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggtree)
library(reshape2)
library(stringr)
library(cowplot)
library(here)

# Check working directory
here()

# Read in tree
tree <- read.tree("./Data/bat_families.nwk")

# Read in data
all_pathogens <- c("Bartonella", "Leptospira", "Mycoplasma", "Rickettsia", "Anaplasma", "Borrelia")

# Empty lists to store data
sampling_sum <- list()
m_sampling_sum <- list()
species_heat_plot <- list()
# Iterate through each pathogen
for(i in 1:length(all_pathogens)){
  # Read in data for each pathogen
  sampling_sum[[i]] <- read_excel("./data/bat_pathogens_review_species_sampling.xlsx", sheet = i)
  # Melt
  m_sampling_sum[[i]] <- sampling_sum[[i]] %>%
    dplyr::select(-Pct_sampled, -Pct_positive) %>%
    melt(id.vars = c("Bat.family", "Plot.order")) %>%
    arrange(Plot.order)
  # Make heatmap image for species sampling
  species_heat_plot[[i]] <-
    ggplot(data = m_sampling_sum[[i]], aes(x = variable, y = reorder(Bat.family, Plot.order))) +
    geom_tile(aes(fill = value), color = "black") +
    geom_text(aes(label = value)) +
    scale_x_discrete(name = "Status") +
    scale_fill_gradientn(name = "Species", colors = c("white", hcl.colors(
      n = 12, palette = "Sunset", rev = TRUE
    ))) +
    theme_cowplot(font_size = 12) +
    theme(
      axis.text.x = element_text(
        size = 10,
        angle = 45,
        hjust = 1
      ),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = 0.25,
      plot.background = element_rect("white")
    )
}

# Empty lists to store data
fam_sum <- list()
m_fam_sum <- list()
fam_heat_plot <- list()
# Iterate through each pathogen
for(i in 1:length(all_pathogens)){
  # Read in data for each pathogen
  fam_sum[[i]] <- read_excel("./data/bat_pathogens_review_family_summary.xlsx", sheet = i)
  # Melt
  m_fam_sum[[i]] <-
    melt(fam_sum, id.vars = c("Bat.family", "Plot.order")) %>%
    arrange(Plot.order)
  # Assign order to value
  m_fam_sum[[i]]$value <- factor(m_fam_sum[[i]]$value, levels = c(4, 3, 2, 1))
  
  # Make heatmap image for family biogeography
  fam_heat_plot[[i]] <-
    ggplot(data = as.data.frame(m_fam_sum[[i]]), aes(
      x = variable,
      y = reorder(Bat.family, Plot.order),
      fill = factor(value)
    )) +
    geom_tile(color = "black") +
    scale_fill_manual(
      name = "Legend",
      values = c(hcl.colors(n = 3, palette = "Sunset"), "#FFFFFF"),
      labels = c(
        "Family present, sampled, pathogen detected",
        "Family present, sampled, pathogen not detected",
        "Family present, not sampled",
        "Family not present"
      )
    ) +
    scale_x_discrete(
      name = "Continent",
      labels = c("N. America", "S. America", "Europe", "Africa", "Asia", "Oceania")
    ) +
    ylab(NULL) +
    theme_cowplot(font_size = 12) +
    theme(
      axis.text.x = element_text(
        size = 10,
        angle = 45,
        hjust = 1
      ),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = 0.5,
      plot.background = element_rect("white")
    )
}

# Iterate through each pathogen
for(i in 1:length(all_pathogens)){
  # Combine images
  panelsAB <- plot_grid(
    fam_heat_plot[[i]],
    species_heat_plot[[i]],
    labels = c("A", "B"),
    label_size = 16,
    ncol = 2,
    rel_widths = c(.55, .45),
    align = "h"
  )
  plot_grid(
    ggtree(tree) + ylim(-1.9, 26.2),
    panelsAB,
    ncol = 2,
    rel_widths = c(.1, .9)
  )
  # Save to file
  ggsave(
    filename = paste0("./results/bat_sampling_summary_", all_pathogens[i], ".pdf"),
    device = "pdf",
    width = 10,
    height = 7,
    units = "in"
  )
  ggsave(
    filename = paste0("./results/bat_sampling_summary_", all_pathogens[i], ".png"),
    device = "png",
    width = 10,
    height = 7,
    dpi = 300,
    units = "in"
  )
  ggsave(
    filename = paste0("./results/bat_sampling_summary_", all_pathogens[i], ".tiff"),
    device = "tiff",
    width = 10,
    height = 7,
    dpi = 300,
    units = "in"
  )
}
