
# PREAMBLE ----------------------------------------------------------------

# Code for meta-analysis of bacterial pathogen prevalence in bats
# R version 4.2.1 "Funny-Looking Kid"
# Works: 2022-11-09

# load libraries
library(here)
library(tidyverse)
library(ape)
library(metafor)

# Define functions
"%ni%" <- Negate("%in%")

# DATA IMPORT -------------------------------------------------------------

# pathogen presence-absence
all_data <- read.csv(here("data", "new_Data_pat_ORIGINAL_2022_updated0310_nocoord.csv"))
# mammalian supertree
mammal_tree <- read.nexus(here("data", "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"))
# mammalian supertree taxonomic data
mammal_taxa <- read.csv(here("data", "taxonomy_mamPhy_5911species.csv"))

# PREPARE DATA ------------------------------------------------------------

# filter pathogen data to valid species
filtered_data <- all_data %>%
  filter(!str_detect(Valid_bat_name_text, coll("cf.")),
         !str_detect(Valid_bat_name_text, coll("sp.")),
         !str_detect(Valid_bat_name_text, coll("spp.")),
         !str_detect(Valid_bat_name_text, coll("idae")),
         !str_detect(Valid_bat_name_text, coll("/")),
         Valid_bat_name_text != "Several",
         for_prevalence_analysis == 1,
         !is.na(Country)) %>%
  drop_na(n_pos, n_tested) %>%
  mutate(Valid_bat_genus = word(Valid_bat_name_text),
         Continent = case_when(Country %in% c("United States",
                                                       "Belize",
                                                       "Costa Rica",
                                                       "Federation of Saint Christopher and Nevis",
                                                       "Grenada",
                                                       "Guatemala",
                                                       "Mexico",
                                                       "Puerto Rico") ~ "North America",
                                        Country %in% c("Argentina",
                                                       "Brazil",
                                                       "Chile",
                                                       "Colombia",
                                                       "French Guiana",
                                                       "Peru") ~ "South America",
                                        Country %in% c("Austria",
                                                       "Czech Republic",
                                                       "Czech Republic/Slovakia",
                                                       "Finland",
                                                       "France",
                                                       "France/Spain",
                                                       "Germany",
                                                       "Hungary",
                                                       "Italy",
                                                       "Multiple",
                                                       "Netherlands",
                                                       "Poland",
                                                       "Romania",
                                                       "Spain",
                                                       "Switzerland",
                                                       "UK") ~ "Europe",
                                        Country %in% c("Tunisia",
                                                       "Comoros",
                                                       "Madagascar",
                                                       "Mauritius",
                                                       "Mayotte",
                                                       "Mozambique",
                                                       "Nigeria",
                                                       "Reunion Island",
                                                       "South Africa",
                                                       "Swaziland",
                                                       "Zambia") ~ "Africa",
                                        Country %in% c("Armenia",
                                                       "China",
                                                       "Georgia",
                                                       "Japan",
                                                       "Russia",
                                                       "Laos",
                                                       "Malaysia",
                                                       "Thailand",
                                                       "Vietnam") ~ "Asia",
                                        Country %in% c("Australia",
                                                       "New Caledonia") ~ "Australia"),
         Ecoregion = case_when(Country %in% c("United States") ~ "Nearctic",
                                        Country %in% c("Austria",
                                                       "Armenia",
                                                       "China",
                                                       "Czech Republic",
                                                       "Czech Republic/Slovakia",
                                                       "Finland",
                                                       "France",
                                                       "France/Spain",
                                                       "Georgia",
                                                       "Germany",
                                                       "Hungary",
                                                       "Italy",
                                                       "Japan",
                                                       "Multiple",
                                                       "Netherlands",
                                                       "Poland",
                                                       "Romania",
                                                       "Russia",
                                                       "Spain",
                                                       "Switzerland",
                                                       "Tunisia",
                                                       "UK") ~ "Palearctic",
                                        Country %in% c("Comoros",
                                                       "Madagascar",
                                                       "Mauritius",
                                                       "Mayotte",
                                                       "Mozambique",
                                                       "Nigeria",
                                                       "Reunion Island",
                                                       "South Africa",
                                                       "Swaziland",
                                                       "Zambia") ~ "Afrotropical",
                                        Country %in% c("Laos",
                                                       "Malaysia",
                                                       "Thailand",
                                                       "Vietnam") ~ "Indomalayan",
                                        Country %in% c("Australia",
                                                       "New Caledonia") ~ "Australasian/Oceanian",
                                        Country %in% c("Argentina",
                                                       "Belize",
                                                       "Brazil",
                                                       "Chile",
                                                       "Colombia",
                                                       "Costa Rica",
                                                       "Federation of Saint Christopher and Nevis",
                                                       "French Guiana",
                                                       "Grenada",
                                                       "Guatemala",
                                                       "Mexico",
                                                       "Peru",
                                                       "Puerto Rico") ~ "Neotropical"))
filt_sum_data <- filtered_data %>%
  group_by(Lit, Pathogen_valid, Bat_family, Valid_bat_name_text, Valid_bat_genus) %>%
  summarize(total_positive = sum(n_pos),
            total_tested = sum(n_tested),
            prevalence = total_positive/total_tested)
filt_sum_data$observation <- factor(1:nrow(filt_sum_data))

# select only bat taxa
bat_taxa <- mammal_taxa[mammal_taxa$ord == "CHIROPTERA",]
bat_taxa <- mammal_taxa[mammal_taxa$ord == "CHIROPTERA",]
bat_taxa$tip <- bat_taxa$Species_Name
bat_taxa$tip <- as.character(bat_taxa$tip)
bat_taxa$fam <- as.character(bat_taxa$fam)

## trim tree to bats
bat_taxa$tiplabel <- as.character(bat_taxa$tiplabel)
bat_tree <- keep.tip(mammal_tree, bat_taxa$tiplabel)

## fix tip labels (remove underscore and replace with space)
bat_tree$tip.label <- sapply(strsplit(bat_tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
bat_taxa$species <- sapply(strsplit(bat_taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## test if all bats are in tree
setdiff(filt_sum_data$Valid_bat_name_text, bat_taxa$species)

#fix bats for dataset being given to reader (e.g. rows with only genus-level info are kept) - not for use in analyses
filt_sum_data$species_for_reader1 <- filt_sum_data$Valid_bat_name_text
filt_sum_data$species_for_reader1 <-plyr::revalue(filt_sum_data$species_for_reader1,
                                c("Diaemus youngii" = "Diaemus youngi",
                                  "Lyroderma lyra" = "Megaderma lyra",
                                  "Eumops nanus" = "Eumops",
                                  "Molossus nigricans" = "Molossus",
                                  "Pteronotus fulvus" = "Pteronotus",
                                  "Pteronotus mesoamericanus" = "Pteronotus",
                                  "Artibeus intermedius" = "Artibeus",
                                  "Gardnerycteris keenani" = "Gardnerycteris",
                                  "Uroderma convexum" = "Uroderma",
                                  "Myotis pilosatibialis" = "Myotis",
                                  "Rhogeessa aenea" = "Rhogeessa aeneus",
                                  "Dermanura gnoma" = "Dermanura gnomus",
                                  "Mops nigeriae" = "Chaerephon nigeriae",
                                  "Dermanura tolteca" = "Dermanura toltecus",
                                  "Hypsugo savii" = "Pipistrellus savii",
                                  "Dermanura cinerea" = "Dermanura cinereus",
                                  "Pteropus vetula" = "Pteropus vetulus",
                                  "Myotis goudotii" = "Myotis goudoti",
                                  "Mops ansorgei" = "Chaerephon ansorgei",
                                  "Mops pumilus" = "Chaerephon pumilus",
                                  "Afronycteris helios" = "Neoromicia helios",
                                  "Afronycteris nana" = "Neoromicia nana",
                                  "Laephotis capensis" = "Neoromicia capensis",
                                  "Vansonia rueppellii" = "Pipistrellus rueppellii",
                                  "Paremballonura tiavato" = "Emballonura tiavato",
                                  "Macronycteris commersonii" = "Hipposideros commersoni",
                                  "Paratriaenops furcula" = "Paratriaenops furculus",
                                  "Mops atsinanana" = "Chaerephon atsinanana",
                                  "Mops leucogaster" = "Mops",
                                  "Laephotis malagasyensis" = "Neoromicia malagasyensis",
                                  "Laephotis matroka" = "Neoromicia matroka",
                                  "Laephotis robertsi" = "Neoromicia robertsi",
                                  "Neoromicia bemainty" = "Hypsugo bemainty",
                                  "Miniopterus orianae bassanii" = "Miniopterus",
                                  "Natalus macrourus" = "Natalus espiritosantensis",
                                  "Gardnerycteris crenulatum" = "Mimon crenulatum",
                                  "Mops pusillus" = "Mops",
                                  "Neoromicia anchieta" = "Hypsugo anchietae",
                                  "Pteropus seychellensis comorensis" = "Pteropus seychellensis",
                                  "Vampyriscus nymphaea" = "Vampyressa nymphaea",
                                  "Mops plicatus" = "Chaerephon plicatus",
                                  "Myotis caucensis" = "Myotis",
                                  "Hipposideros gentilis" = "Hipposideros",
                                  "Hipposideros kunzi" = "Hipposideros",
                                  "Rhinolophus refulgens" = "Rhinolophus",
                                  "Myotis sibiricus" = "Myotis",
                                  "Penthetor lucasii" = "Penthetor lucasi"))

# make an alternative dataset where species are replaces with their most closely related species
filt_sum_data$species_for_reader2 <- filt_sum_data$Valid_bat_name_text
filt_sum_data$species_for_reader2 <-plyr::revalue(filt_sum_data$species_for_reader2,
                                                  c("Diaemus youngii" = "Diaemus youngi",
                                                    "Lyroderma lyra" = "Megaderma lyra",
                                                    "Eumops nanus" = "Eumops bonariensis",
                                                    "Molossus nigricans" = "Molossus rufus",
                                                    "Pteronotus fulvus" = "Pteronotus davyi",
                                                    "Pteronotus mesoamericanus" = "Pteronotus parnellii",
                                                    "Artibeus intermedius" = "Artibeus lituratus",
                                                    "Gardnerycteris keenani" = "Mimon crenulatum",
                                                    "Uroderma convexum" = "Uroderma bilobatum",
                                                    "Myotis pilosatibialis" = "Myotis keaysi",
                                                    "Rhogeessa aenea" = "Rhogeessa aeneus",
                                                    "Dermanura gnoma" = "Dermanura gnomus",
                                                    "Mops nigeriae" = "Chaerephon nigeriae",
                                                    "Dermanura tolteca" = "Dermanura toltecus",
                                                    "Hypsugo savii" = "Pipistrellus savii",
                                                    "Dermanura cinerea" = "Dermanura cinereus",
                                                    "Pteropus vetula" = "Pteropus vetulus",
                                                    "Myotis goudotii" = "Myotis goudoti",
                                                    "Mops ansorgei" = "Chaerephon ansorgei",
                                                    "Mops pumilus" = "Chaerephon pumilus",
                                                    "Afronycteris helios" = "Neoromicia helios",
                                                    "Afronycteris nana" = "Neoromicia nana",
                                                    "Laephotis capensis" = "Neoromicia capensis",
                                                    "Vansonia rueppellii" = "Pipistrellus rueppellii",
                                                    "Paremballonura tiavato" = "Emballonura tiavato",
                                                    "Macronycteris commersonii" = "Hipposideros commersoni",
                                                    "Paratriaenops furcula" = "Paratriaenops furculus",
                                                    "Mops atsinanana" = "Chaerephon atsinanana",
                                                    "Mops leucogaster" = "Chaerephon pumilus",
                                                    "Laephotis malagasyensis" = "Neoromicia malagasyensis",
                                                    "Laephotis matroka" = "Neoromicia matroka",
                                                    "Laephotis robertsi" = "Neoromicia robertsi",
                                                    "Neoromicia bemainty" = "Hypsugo bemainty",
                                                    "Miniopterus orianae bassanii" = "Miniopterus schreibersii",
                                                    "Natalus macrourus" = "Natalus espiritosantensis",
                                                    "Gardnerycteris crenulatum" = "Mimon crenulatum",
                                                    "Mops pusillus" = "Chaerephon pumilus",
                                                    "Neoromicia anchieta" = "Hypsugo anchietae",
                                                    "Pteropus seychellensis comorensis" = "Pteropus seychellensis",
                                                    "Vampyriscus nymphaea" = "Vampyressa nymphaea",
                                                    "Mops plicatus" = "Chaerephon plicatus",
                                                    "Myotis caucensis" = "Myotis nigricans",
                                                    "Hipposideros gentilis" = "Hipposideros pomona",
                                                    "Hipposideros kunzi" = "Hipposideros bicolor",
                                                    "Rhinolophus refulgens" = "Rhinolophus lepidus",
                                                    "Myotis sibiricus" = "Myotis brandtii",
                                                    "Penthetor lucasii" = "Penthetor lucasi"))

# check
setdiff(filt_sum_data$species_for_reader1, bat_tree$tip.label)
setdiff(filt_sum_data$species_for_reader1, bat_taxa$species)

setdiff(filt_sum_data$species_for_reader2, bat_tree$tip.label)
setdiff(filt_sum_data$species_for_reader2, bat_taxa$species)

# remove genus-only/drop rows for analyses (trim dataset to only phylogeny)
dataset1 <- filt_sum_data[filt_sum_data$species_for_reader1 %in% bat_tree$tip.label,] %>%
  select(-species_for_reader2) %>%
  rename(species_for_reader = species_for_reader1)
# alternative dataset with replacement species
dataset2 <- filt_sum_data[filt_sum_data$species_for_reader2 %in% bat_tree$tip.label,] %>%
  select(-species_for_reader1) %>%
  rename(species_for_reader = species_for_reader2)

# FILTER DATA AND RUN MODELS ----------------------------------------------

# list of pathogens
pathogen_list <- c("Bartonella", "Leptospira", "Mycoplasma",
                   "Rickettsia", "Anaplasma", "Borrelia",
                   "Coxiella")
# function for I2 for rma.mv
i2=function(model){
  
  # metafor site code for I2
  W=diag(1/model$vi)
  X=model.matrix(model)
  P=W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2=100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2=round(I2,2)
  
  # summarize by each variance component
  allI2=100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  allI2=round(allI2,3)
  return(list(I2=I2,allI2=allI2))
}

# function to filter data to pathogen and run GLMM
run_model <- function(dataset, pathogen, minimum_tested, model_formula){
  # filter data to pathogen of interest and to the families with sufficient data
  filtered_dataset <- dataset %>%
    filter(Pathogen_valid == pathogen)
  filtered_dataset_lowfamily <- filtered_dataset %>%
    group_by(Bat_family) %>%
    summarize(total_tested2 = sum(total_tested)) %>%
    filter(total_tested2 < minimum_tested)
  filtered_dataset <- filtered_dataset %>%
    filter(Bat_family %ni% unique(filtered_dataset_lowfamily$Bat_family))
  
  # calculate pft in escalc for yi and vi
  filtered_dataset <- data.frame(filtered_dataset,
                                    escalc(xi = filtered_dataset$total_positive,
                                           ni = filtered_dataset$total_tested,
                                           measure="PFT"))
  # back transform
  filtered_dataset$backtrans <- transf.ipft(filtered_dataset$yi, filtered_dataset$total_tested)
  
  # species and phylo effect
  filtered_dataset$phylo <- filtered_dataset$species_for_reader
  filtered_dataset$species <- filtered_dataset$phylo
  
  # fix tips in tree
  stree_filtered_dataset <- keep.tip(bat_tree, as.character(unique(filtered_dataset$species)))
  
  # convert tree to correlation matrix
  cmatrix_filtered_dataset <- vcv.phylo(stree_filtered_dataset, cor=T)
  
  # run model
  model_run <- rma.mv(yi = yi,
                              V = vi,
                              random = list(~1|Lit, ~1|Valid_bat_name_text, ~1|phylo),
                              R = list(phylo = cmatrix_filtered_dataset),
                              method = "REML",
                              mods = model_formula,
                              data = filtered_dataset,
                              control = list(optimizer = "optim", optmethod = "BFGS"))
  
  model_run_alt <- rma.mv(yi = yi,
                          V = vi,
                          random = list(~1|Lit, ~1|Valid_bat_name_text),
                          method = "REML",
                          mods = model_formula,
                          data = filtered_dataset,
                          control = list(optimizer = "optim", optmethod = "BFGS"))
  
  # output data
  return(list(data = filtered_dataset, model = model_run, model_alt = model_run_alt))
}

null_model_output1 <- list(NULL)
null_model_output2 <- list(NULL)
model_output1 <- list(NULL)
model_output2 <- list(NULL)
prev_model_output1 <- list(NULL)
prev_model_output2 <- list(NULL)
formula_option1 <- as.formula("~Bat_family") # with intercept for test of significance for bat families
formula_option2 <- as.formula("~Bat_family - 1") # without intercept for easy output of estimated prevalence with confidence intervals
for(i in 1:length(pathogen_list)){
  null_model_output1[[i]] <- run_model(dataset1, pathogen_list[i], 1, as.formula("~1"))
  null_model_output2[[i]] <- run_model(dataset2, pathogen_list[i], 1, as.formula("~1"))
  model_output1[[i]] <- run_model(dataset1, pathogen_list[i], 1, formula_option1)
  model_output2[[i]] <- run_model(dataset2, pathogen_list[i], 1, formula_option1)
  prev_model_output1[[i]] <- run_model(dataset1, pathogen_list[i], 1, formula_option2)
  prev_model_output2[[i]] <- run_model(dataset2, pathogen_list[i], 1, formula_option2)
}

# FORMATTING MODEL OUTPUT -------------------------------------------------

# ~ phylogenetic meta-analysis --------------------------------------------

# inspect model fits, variance, and heterogeneity
# model fits with where 13 species without match in phylogeny were removed
# Bartonella
summary(null_model_output1[[1]]$model)
i2(null_model_output1[[1]]$model)
(deviance(null_model_output1[[1]]$model) - deviance(model_output1[[1]]$model)) / deviance(null_model_output1[[1]]$model)
summary(model_output1[[1]]$model)
var(model_output1[[1]]$data$prevalence)
# Leptospira
summary(null_model_output1[[2]]$model)
i2(null_model_output1[[2]]$model)
(deviance(null_model_output1[[2]]$model) - deviance(model_output1[[2]]$model)) / deviance(null_model_output1[[2]]$model)
summary(model_output1[[2]]$model)
var(model_output1[[2]]$data$prevalence)
# Mycoplasma
summary(null_model_output1[[3]]$model)
i2(null_model_output1[[3]]$model)
(deviance(null_model_output1[[3]]$model) - deviance(model_output1[[3]]$model)) / deviance(null_model_output1[[3]]$model)
summary(model_output1[[3]]$model)
var(model_output1[[3]]$data$prevalence)
# Rickettsia
summary(null_model_output1[[4]]$model)
i2(null_model_output1[[4]]$model)
(deviance(null_model_output1[[4]]$model) - deviance(model_output1[[4]]$model)) / deviance(null_model_output1[[4]]$model)
summary(model_output1[[4]]$model)
var(model_output1[[4]]$data$prevalence)
# Anaplasma
summary(null_model_output1[[5]]$model)
i2(null_model_output1[[5]]$model)
(deviance(null_model_output1[[5]]$model) - deviance(model_output1[[5]]$model)) / deviance(null_model_output1[[5]]$model)
summary(model_output1[[5]]$model)
var(model_output1[[5]]$data$prevalence)
# Borrelia
summary(null_model_output1[[6]]$model)
i2(null_model_output1[[6]]$model)
(deviance(null_model_output1[[6]]$model) - deviance(model_output1[[6]]$model)) / deviance(null_model_output1[[6]]$model)
summary(model_output1[[6]]$model)
var(model_output1[[6]]$data$prevalence)
# Coxiella
summary(null_model_output1[[7]]$model)
i2(null_model_output1[[7]]$model)
(deviance(null_model_output1[[7]]$model) - deviance(model_output1[[7]]$model)) / deviance(null_model_output1[[7]]$model)
summary(model_output1[[7]]$model)
var(model_output1[[7]]$data$prevalence)

# model fits with replacements for 13 species
# Bartonella
summary(null_model_output2[[1]]$model)
i2(null_model_output2[[1]]$model)
(deviance(null_model_output2[[1]]$model) - deviance(model_output2[[1]]$model)) / deviance(null_model_output2[[1]]$model)
summary(model_output2[[1]]$model)
var(model_output2[[1]]$data$prevalence)
# Leptospira
summary(null_model_output2[[2]]$model)
i2(null_model_output2[[2]]$model)
(deviance(null_model_output2[[2]]$model) - deviance(model_output2[[2]]$model)) / deviance(null_model_output2[[2]]$model)
summary(model_output2[[2]]$model)
var(model_output2[[2]]$data$prevalence)
# Mycoplasma
summary(null_model_output2[[3]]$model)
i2(null_model_output2[[3]]$model)
(deviance(null_model_output2[[3]]$model) - deviance(model_output2[[3]]$model)) / deviance(null_model_output2[[3]]$model)
summary(model_output2[[3]]$model)
var(model_output2[[3]]$data$prevalence)
# Rickettsia
summary(null_model_output2[[4]]$model)
i2(null_model_output2[[4]]$model)
(deviance(null_model_output2[[4]]$model) - deviance(model_output2[[4]]$model)) / deviance(null_model_output2[[4]]$model)
summary(model_output2[[4]]$model)
var(model_output2[[4]]$data$prevalence)
# Anaplasma
summary(null_model_output2[[5]]$model)
i2(null_model_output2[[5]]$model)
(deviance(null_model_output2[[5]]$model) - deviance(model_output2[[5]]$model)) / deviance(null_model_output2[[5]]$model)
summary(model_output2[[5]]$model)
var(model_output2[[5]]$data$prevalence)
# Borrelia
summary(null_model_output2[[6]]$model)
i2(null_model_output2[[6]]$model)
(deviance(null_model_output2[[6]]$model) - deviance(model_output2[[6]]$model)) / deviance(null_model_output2[[6]]$model)
summary(model_output2[[6]]$model)
var(model_output2[[6]]$data$prevalence)
# Coxiella
summary(null_model_output2[[7]]$model)
i2(null_model_output2[[7]]$model)
(deviance(null_model_output2[[7]]$model) - deviance(model_output2[[7]]$model)) / deviance(null_model_output2[[7]]$model)
summary(model_output2[[7]]$model)
var(model_output2[[7]]$data$prevalence)

# loop through pathogens to calculate prevalence and meta-analysis fitted prevalence
for(i in 1:length(pathogen_list)){
  prevalence <- all_data %>%
    filter(Pathogen_valid == pathogen_list[i],
           for_prevalence_analysis == 1) %>%
    drop_na(n_pos, n_tested) %>%
    group_by(Bat_family) %>%
    summarize(Grand_total_positive = sum(n_pos),
              Grand_total_tested = sum(n_tested)) %>%
    mutate(binom::binom.exact(x = Grand_total_positive, n = Grand_total_tested),
           mean = round(mean, 3)*100,
           lower = round(lower, 3)*100,
           upper = round(upper, 3)*100)
  meta_prevalence <- prev_model_output2[[i]]$data %>%
    group_by(Bat_family) %>%
    summarize(Sub_total_positive = sum(total_positive),
              Sub_total_tested = sum(total_tested))
  meta_prevalence <- data.frame(meta_prevalence,
                                           beta = prev_model_output2[[i]]$model$beta,
                                           ci.lb = prev_model_output2[[i]]$model$ci.lb,
                                           ci.ub = prev_model_output2[[i]]$model$ci.ub) %>%
    mutate(fitted.prev = round(transf.ipft(beta, Sub_total_tested), 3)*100,
           fitted.lb = round(transf.ipft(ci.lb, Sub_total_tested), 3)*100,
           fitted.ub = round(transf.ipft(ci.ub, Sub_total_tested), 3)*100)
  output <- left_join(prevalence, meta_prevalence) %>%
    mutate(concat_prev = paste0(Grand_total_positive, "/", Grand_total_tested, " (", mean, " [", lower, ", ", upper, "])"),
           concat_meta = paste0(Sub_total_positive, "/", Sub_total_tested, " (", fitted.prev, " [", fitted.lb, ", ", fitted.ub, "])"))
  write.csv(output, here("results", paste0(pathogen_list[i], "_fitted_prevalence.csv")), row.names = FALSE)
}

# check prevalence
check_prevalence <- all_data %>%
  filter(for_prevalence_analysis == 1) %>%
  drop_na(n_pos, n_tested) %>%
  group_by(Pathogen_valid, Bat_family) %>%
  summarize(Grand_total_positive = sum(n_pos),
            Grand_total_tested = sum(n_tested)) %>%
  mutate(binom::binom.exact(x = Grand_total_positive, n = Grand_total_tested),
         mean = round(mean, 3)*100,
         lower = round(lower, 3)*100,
         upper = round(upper, 3)*100) %>%
  mutate(concat_prev = paste0(Grand_total_positive, "/", Grand_total_tested, " (", mean, " [", lower, ", ", upper, "])"))
View(check_prevalence)

# ~ regular meta-analysis (no phylogenetic random effect) ------------------

# inspect model fits, variance, and heterogeneity
# model fits with where 13 species without match in phylogeny were removed
# Bartonella
summary(null_model_output1[[1]]$model_alt)
i2(null_model_output1[[1]]$model_alt)
(deviance(null_model_output1[[1]]$model_alt) - deviance(model_output1[[1]]$model_alt)) / deviance(null_model_output1[[1]]$model_alt)
summary(model_output1[[1]]$model_alt)
var(model_output1[[1]]$data$prevalence)
# Leptospira
summary(null_model_output1[[2]]$model_alt)
i2(null_model_output1[[2]]$model_alt)
(deviance(null_model_output1[[2]]$model_alt) - deviance(model_output1[[2]]$model_alt)) / deviance(null_model_output1[[2]]$model_alt)
summary(model_output1[[2]]$model_alt)
var(model_output1[[2]]$data$prevalence)
# Mycoplasma
summary(null_model_output1[[3]]$model_alt)
i2(null_model_output1[[3]]$model_alt)
(deviance(null_model_output1[[3]]$model_alt) - deviance(model_output1[[3]]$model_alt)) / deviance(null_model_output1[[3]]$model_alt)
summary(model_output1[[3]]$model_alt)
var(model_output1[[3]]$data$prevalence)
# Rickettsia
summary(null_model_output1[[4]]$model_alt)
i2(null_model_output1[[4]]$model_alt)
(deviance(null_model_output1[[4]]$model_alt) - deviance(model_output1[[4]]$model_alt)) / deviance(null_model_output1[[4]]$model_alt)
summary(model_output1[[4]]$model_alt)
var(model_output1[[4]]$data$prevalence)
# Anaplasma
summary(null_model_output1[[5]]$model_alt)
i2(null_model_output1[[5]]$model_alt)
(deviance(null_model_output1[[5]]$model_alt) - deviance(model_output1[[5]]$model_alt)) / deviance(null_model_output1[[5]]$model_alt)
summary(model_output1[[5]]$model_alt)
var(model_output1[[5]]$data$prevalence)
# Borrelia
summary(null_model_output1[[6]]$model_alt)
i2(null_model_output1[[6]]$model_alt)
(deviance(null_model_output1[[6]]$model_alt) - deviance(model_output1[[6]]$model_alt)) / deviance(null_model_output1[[6]]$model_alt)
summary(model_output1[[6]]$model_alt)
var(model_output1[[6]]$data$prevalence)
# Coxiella
summary(null_model_output1[[7]]$model_alt)
i2(null_model_output1[[7]]$model_alt)
(deviance(null_model_output1[[7]]$model_alt) - deviance(model_output1[[7]]$model_alt)) / deviance(null_model_output1[[7]]$model_alt)
summary(model_output1[[7]]$model_alt)
var(model_output1[[7]]$data$prevalence)

# model fits with replacements for 13 species
# Bartonella
summary(null_model_output2[[1]]$model_alt)
i2(null_model_output2[[1]]$model_alt)
(deviance(null_model_output2[[1]]$model_alt) - deviance(model_output2[[1]]$model_alt)) / deviance(null_model_output2[[1]]$model_alt)
summary(model_output2[[1]]$model_alt)
var(model_output2[[1]]$data$prevalence)
# Leptospira
summary(null_model_output2[[2]]$model_alt)
i2(null_model_output2[[2]]$model_alt)
(deviance(null_model_output2[[2]]$model_alt) - deviance(model_output2[[2]]$model_alt)) / deviance(null_model_output2[[2]]$model_alt)
summary(model_output2[[2]]$model_alt)
var(model_output2[[2]]$data$prevalence)
# Mycoplasma
summary(null_model_output2[[3]]$model_alt)
i2(null_model_output2[[3]]$model_alt)
(deviance(null_model_output2[[3]]$model_alt) - deviance(model_output2[[3]]$model_alt)) / deviance(null_model_output2[[3]]$model_alt)
summary(model_output2[[3]]$model_alt)
var(model_output2[[3]]$data$prevalence)
# Rickettsia
summary(null_model_output2[[4]]$model_alt)
i2(null_model_output2[[4]]$model_alt)
(deviance(null_model_output2[[4]]$model_alt) - deviance(model_output2[[4]]$model_alt)) / deviance(null_model_output2[[4]]$model_alt)
summary(model_output2[[4]]$model_alt)
var(model_output2[[4]]$data$prevalence)
# Anaplasma
summary(null_model_output2[[5]]$model_alt)
i2(null_model_output2[[5]]$model_alt)
(deviance(null_model_output2[[5]]$model_alt) - deviance(model_output2[[5]]$model_alt)) / deviance(null_model_output2[[5]]$model_alt)
summary(model_output2[[5]]$model_alt)
var(model_output2[[5]]$data$prevalence)
# Borrelia
summary(null_model_output2[[6]]$model_alt)
i2(null_model_output2[[6]]$model_alt)
(deviance(null_model_output2[[6]]$model_alt) - deviance(model_output2[[6]]$model_alt)) / deviance(null_model_output2[[6]]$model_alt)
summary(model_output2[[6]]$model_alt)
var(model_output2[[6]]$data$prevalence)
# Coxiella
summary(null_model_output2[[7]]$model_alt)
i2(null_model_output2[[7]]$model_alt)
(deviance(null_model_output2[[7]]$model_alt) - deviance(model_output2[[7]]$model_alt)) / deviance(null_model_output2[[7]]$model_alt)
summary(model_output2[[7]]$model_alt)
var(model_output2[[7]]$data$prevalence)

# loop through pathogens to calculate prevalence and meta-analysis fitted prevalence
for(i in 1:length(pathogen_list)){
  prevalence_alt <- all_data %>%
    filter(Pathogen_valid == pathogen_list[i],
           for_prevalence_analysis == 1) %>%
    drop_na(n_pos, n_tested) %>%
    group_by(Bat_family) %>%
    summarize(Grand_total_positive = sum(n_pos),
              Grand_total_tested = sum(n_tested)) %>%
    mutate(binom::binom.exact(x = Grand_total_positive, n = Grand_total_tested),
           mean = round(mean, 3)*100,
           lower = round(lower, 3)*100,
           upper = round(upper, 3)*100)
  meta_prevalence <- prev_model_output2[[i]]$data %>%
    group_by(Bat_family) %>%
    summarize(Sub_total_positive = sum(total_positive),
              Sub_total_tested = sum(total_tested))
  meta_prevalence <- data.frame(meta_prevalence,
                                beta = prev_model_output2[[i]]$model_alt$beta,
                                ci.lb = prev_model_output2[[i]]$model_alt$ci.lb,
                                ci.ub = prev_model_output2[[i]]$model_alt$ci.ub) %>%
    mutate(fitted.prev = round(transf.ipft(beta, Sub_total_tested), 3)*100,
           fitted.lb = round(transf.ipft(ci.lb, Sub_total_tested), 3)*100,
           fitted.ub = round(transf.ipft(ci.ub, Sub_total_tested), 3)*100)
  output <- left_join(prevalence, meta_prevalence) %>%
    mutate(concat_prev = paste0(Grand_total_positive, "/", Grand_total_tested, " (", mean, " [", lower, ", ", upper, "])"),
           concat_meta = paste0(Sub_total_positive, "/", Sub_total_tested, " (", fitted.prev, " [", fitted.lb, ", ", fitted.ub, "])"))
  write.csv(output, here("results", paste0(pathogen_list[i], "_fitted_prevalence_alt.csv")), row.names = FALSE)
}

# check prevalence
check_prevalence_alt <- all_data %>%
  filter(for_prevalence_analysis == 1) %>%
  drop_na(n_pos, n_tested) %>%
  group_by(Pathogen_valid, Bat_family) %>%
  summarize(Grand_total_positive = sum(n_pos),
            Grand_total_tested = sum(n_tested)) %>%
  mutate(binom::binom.exact(x = Grand_total_positive, n = Grand_total_tested),
         mean = round(mean, 3)*100,
         lower = round(lower, 3)*100,
         upper = round(upper, 3)*100) %>%
  mutate(concat_prev = paste0(Grand_total_positive, "/", Grand_total_tested, " (", mean, " [", lower, ", ", upper, "])"))
View(check_prevalence_alt)
