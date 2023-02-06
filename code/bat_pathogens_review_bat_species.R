##########################
# Code for mapping overlap between bat species
# R version 4.2.1 "Funny-Looking Kid"
# Works: 2022-10-07
##########################

# Filter to only the order Chiroptera
myspecies.distr <- mammterr[mammterr$order_ == "CHIROPTERA",]

##################################
### Species distribution plots ###
##################################

# Sample points across the species distributions on a regular grid
spp.points <-
  spsample(myspecies.distr, n = 1000000, type = "regular")
# spp.points <-
#   spsample(myspecies.distr, n = 1000, type = "regular")

# Summarize the species distributions that overlap the chosen points
spp.pointsInPolygons <-
  sp::over(x = spp.points, y = myspecies.distr, returnList = TRUE)

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

# Plot map of species range overlap for all coronavirus hosts
plotA <- ggplot() +
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
  ggtitle("Described and mapped bat species") +
  theme(
    plot.title = element_text(size = 24, face = 2),
    plot.background = element_rect(fill = "white", color = NA)
  )

###################
### End of code ###
###################
