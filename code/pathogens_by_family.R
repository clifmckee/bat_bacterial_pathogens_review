tested_species <- bat_spp_sorted %>%
  filter(!str_detect(Valid_bat_name_text, coll("cf.")),
         !str_detect(Valid_bat_name_text, coll("sp.")),
         !str_detect(Valid_bat_name_text, coll("spp.")),
         !str_detect(Valid_bat_name_text, coll("idae")),
         !str_detect(Valid_bat_name_text, coll("/")),
         Valid_bat_name_text != "Several") %>%
  group_by(Bat_family, Pathogen_valid) %>%
  summarize(
    species_sampled = n_distinct(Valid_bat_name_text)
  )

positive_species <- bat_spp_sorted %>%
  filter(!str_detect(Valid_bat_name_text, coll("cf.")),
         !str_detect(Valid_bat_name_text, coll("sp.")),
         !str_detect(Valid_bat_name_text, coll("spp.")),
         !str_detect(Valid_bat_name_text, coll("idae")),
         !str_detect(Valid_bat_name_text, coll("/")),
         Valid_bat_name_text != "Several",
         n_pos != 0) %>%
  group_by(Bat_family, Pathogen_valid) %>%
  summarize(
    species_tested = n_distinct(Valid_bat_name_text)
  )

merged <- full_join(tested_species, positive_species)
