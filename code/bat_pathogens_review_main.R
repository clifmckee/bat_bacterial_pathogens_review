##########################
# Code for mapping overlap between bat species sampled for bacterial pathogens
# R version 4.2.1 "Funny-Looking Kid"
# Works: 2022-10-07
##########################

# Read in bacterial pathogen and bat species data
bat_spp <-
  read.csv("./data/new_Data_pat_ORIGINAL_2022_updated0310_nocoord.csv")
all_pathogens <-
  c(
    "Bartonella",
    "Leptospira",
    "Mycoplasma",
    "Rickettsia",
    "Anaplasma",
    "Borrelia"
  )

# Sort data
bat_spp_sorted <- bat_spp %>%
  arrange(Bat_family, Valid_bat_name_text, Pathogen)
write.csv(bat_spp_sorted, "./data/sampled_spp_sorted.csv")

# Summarize all data
filtered_bat_spp_sorted <- bat_spp_sorted %>%
  filter(!str_detect(Valid_bat_name_text, coll("cf.")),
         !str_detect(Valid_bat_name_text, coll("sp.")),
         !str_detect(Valid_bat_name_text, coll("spp.")),
         !str_detect(Valid_bat_name_text, coll("idae")),
         !str_detect(Valid_bat_name_text, coll("/")),
         Valid_bat_name_text != "Several") %>%
  mutate(Valid_bat_genus = word(Valid_bat_name_text))

filtered_bat_spp_sorted %>%
  summarize(
    unique.species = n_distinct(Valid_bat_name_text),
    unique.genera = n_distinct(Valid_bat_genus),
    unique.families = n_distinct(Bat_family),
    unique.refs = n_distinct(Lit)
  )

# Summarize all data by family
filtered_bat_spp_sorted %>%
  group_by(Bat_family) %>%
  summarize(unique.species = n_distinct(Valid_bat_name_text))

# Summarize data, separated by pathogen genera
filtered_bat_spp_sorted %>%
  group_by(Pathogen_valid) %>%
  summarize(
    unique.species = n_distinct(Valid_bat_name_text),
    unique.families = n_distinct(Bat_family),
    unique.refs = n_distinct(Lit)
  )

# Correct some names as they appear in the IUCN records
bat_spp[bat_spp$Valid_bat_name_text == "Mops ansorgei", ]$Valid_bat_name_text <-
  "Chaerephon ansorgei"
bat_spp[bat_spp$Valid_bat_name_text == "Mops atsinanana", ]$Valid_bat_name_text <-
  "Chaerephon atsinanana"
bat_spp[bat_spp$Valid_bat_name_text == "Mops nigeriae", ]$Valid_bat_name_text <-
  "Chaerephon nigeriae"
bat_spp[bat_spp$Valid_bat_name_text == "Mops plicatus", ]$Valid_bat_name_text <-
  "Chaerephon plicatus"
bat_spp[bat_spp$Valid_bat_name_text == "Mops pumilus", ]$Valid_bat_name_text <-
  "Chaerephon pumilus"
bat_spp[bat_spp$Valid_bat_name_text == "Diaemus youngii", ]$Valid_bat_name_text <-
  "Diaemus youngi"
bat_spp[bat_spp$Valid_bat_name_text == "Macronycteris commersonii", ]$Valid_bat_name_text <-
  "Macronycteris commersoni"
bat_spp[bat_spp$Valid_bat_name_text == "Neoromicia anchieta", ]$Valid_bat_name_text <-
  "Pipistrellus anchietae"
bat_spp[bat_spp$Valid_bat_name_text == "Neoromicia bemainty", ]$Valid_bat_name_text <-
  "Hypsugo bemainty"
bat_spp[bat_spp$Valid_bat_name_text == "Myonycteris angolensis", ]$Valid_bat_name_text <-
  "Lissonycteris angolensis"
bat_spp[bat_spp$Valid_bat_name_text == "Hsunycteris thomasi", ]$Valid_bat_name_text <-
  "Lonchophylla thomasi"
bat_spp[bat_spp$Valid_bat_name_text == "Natalus macrourus", ]$Valid_bat_name_text <-
  "Natalus espiritosantensis"
bat_spp[bat_spp$Valid_bat_name_text == "Laephotis capensis", ]$Valid_bat_name_text <-
  "Neoromicia capensis"
bat_spp[bat_spp$Valid_bat_name_text == "Afronycteris helios", ]$Valid_bat_name_text <-
  "Neoromicia helios"
bat_spp[bat_spp$Valid_bat_name_text == "Laephotis malagasyensis", ]$Valid_bat_name_text <-
  "Neoromicia malagasyensis"
bat_spp[bat_spp$Valid_bat_name_text == "Laephotis matroka", ]$Valid_bat_name_text <-
  "Neoromicia matroka"
bat_spp[bat_spp$Valid_bat_name_text == "Afronycteris nana", ]$Valid_bat_name_text <-
  "Neoromicia nana"
bat_spp[bat_spp$Valid_bat_name_text == "Laephotis robertsi", ]$Valid_bat_name_text <-
  "Neoromicia robertsi"
bat_spp[bat_spp$Valid_bat_name_text == "Paratriaenops furcula", ]$Valid_bat_name_text <-
  "Paratriaenops furculus"
bat_spp[bat_spp$Valid_bat_name_text == "Penthetor lucasii", ]$Valid_bat_name_text <-
  "Penthetor lucasi"
bat_spp[bat_spp$Valid_bat_name_text == "Pteropus vetula", ]$Valid_bat_name_text <-
  "Pteropus vetulus"
bat_spp[bat_spp$Valid_bat_name_text == "Rhogeessa aenea", ]$Valid_bat_name_text <-
  "Rhogeessa aeneus"

# Find bat species names that do not match IUCN records
bat_spp_missing <- bat_spp %>%
  mutate(IUCN = Valid_bat_name_text %in% mammterr$binomial) %>%
  filter(IUCN == "FALSE")

# Replace some NA values in the n_pos column, create a tested column, and convert n_pos to binary
bat_spp <- bat_spp %>%
  mutate(n_pos = replace_na(n_pos, 1),
         tested = 1,
         n_pos_binary = as.numeric(n_pos > 0))

pathogen_map <- function(pathogen, result, size) {
  # Filter out bad bat species and choose pathogen genera
  unique_bat_spp <- bat_spp %>%
    filter(Valid_bat_name_text %ni% bat_spp_missing$Valid_bat_name_text,
           Pathogen_valid %in% pathogen,
           n_pos_binary %in% result) %>%
    distinct(Valid_bat_name_text) %>%
    arrange(Valid_bat_name_text)
  
  # Create a list of species names to filter mammterr
  allnames = as.character(unique_bat_spp$Valid_bat_name_text)
  
  # Filter mammterr down to just binomial names
  all.binomial = mammterr$binomial
  
  # Check to see if allnames are in the binomial names
  unique.allnames <- unique(allnames)
  
  # Which rows of the data are the species I care about
  keep = list()
  for (i in 1:length(unique.allnames)) {
    keep[[i]] = which(all.binomial == unique.allnames[i])
  }
  x = keep[[1]]
  for (i in 2:length(unique.allnames)) {
    x = c(x, keep[[i]])
  }
  keep.species = data.frame(x = x, species = mammterr[x,]$binomial)
  myspecies.distr <- mammterr[x,]
  
  # Sample points across the species distributions on a regular grid
  spp.points <-
    spsample(myspecies.distr, type = "regular", n = size)
  
  # Summarize the species distributions that overlap the chosen points
  spp.pointsInPolygons <-
    sp::over(x = spp.points,
             y = myspecies.distr,
             returnList = TRUE)
  
  # Count the number of intersections
  spp.counting <-
    lapply(
      spp.pointsInPolygons,
      FUN = function(x)
        nrow(x)
    )
  spp.over.df <-
    data.frame(
      "point" = rownames(t(do.call(
        "cbind", spp.counting
      ))),
      "count" = t(do.call("cbind", spp.counting)),
      "polygon" = paste(spp.pointsInPolygons)
    )
  
  # Summarize counts in a data frame
  spp.points.df <- as.data.frame(spp.points)
  spp.points.df$count <- spp.over.df$count
  
  # Data for world map
  world_map <- map_data("world")
  
  # Combine world map data and species range data into one object
  sf_data <- list(world = world_map, overlap = spp.points.df)
  
  # Plot map of species range overlap
  out <- ggplot() +
    geom_polygon(
      data = sf_data$world,
      aes(x = long, y = lat, group = group),
      color = NA,
      fill = "grey"
    ) +
    geom_tile(data = sf_data$overlap, aes(x = x1, y = x2, fill = count)) +
    ylim(-56, 84) +
    scale_fill_viridis(option = "D",
                       name = "Species",
                       breaks = c(1, ceiling(
                         max(sf_data$overlap$count) * c(0.25, 0.5, 0.75, 1)
                       ))) +
    theme_map(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  return(out)
}

# Plot separate maps of species range overlap for key bacterial pathogen genera hosts
tested_list <- list()
positive_list <- list()
size_list <-
  c(100000, 100000, 100000, 100000, 100000, 100000)
# size_list <-
#   c(1000, 1000, 1000, 1000, 1000, 1000)
for (i in 1:length(all_pathogens)) {
  tested_list[[i]] <-
    pathogen_map(pathogen = all_pathogens[i], result = c(0, 1), size = size_list[i])
  positive_list[[i]] <-
    pathogen_map(pathogen = all_pathogens[i], result = 1, size = size_list[i])
}

plotB <-
  plot_grid(
    # Bartonella, tested
    tested_list[[1]] +
      ggtitle(paste("Species tested for", all_pathogens[1])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Leptospira, tested
    tested_list[[2]] +
      ggtitle(paste("Species tested for", all_pathogens[2])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Mycoplasma, tested
    tested_list[[3]] +
      ggtitle(paste("Species tested for", all_pathogens[3])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Bartonella, positive
    positive_list[[1]] +
      ggtitle(paste("Species positive for", all_pathogens[1])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Leptospira, positive
    positive_list[[2]] +
      ggtitle(paste("Species positive for", all_pathogens[2])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Mycoplasma, positive
    positive_list[[3]] +
      ggtitle(paste("Species positive for", all_pathogens[3])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Rickettsia, tested
    tested_list[[4]] +
      ggtitle(paste("Species tested for", all_pathogens[4])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Anaplasma, tested
    tested_list[[5]] +
      ggtitle(paste("Species tested for", all_pathogens[5])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Borrelia, tested
    tested_list[[6]] +
      ggtitle(paste("Species tested for", all_pathogens[6])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Rickettsia, positive
    positive_list[[4]] +
      ggtitle(paste("Species positive for", all_pathogens[4])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Anaplasma, positive
    positive_list[[5]] +
      ggtitle(paste("Species positive for", all_pathogens[5])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    # Borrelia, positive
    positive_list[[6]] +
      ggtitle(paste("Species positive for", all_pathogens[6])) +
      theme(plot.title = element_text(size = 24, face = 2)),
    nrow = 4,
    ncol = 3
  )

###################
### End of code ###
###################
