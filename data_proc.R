library(tidyverse)
library(lubridate)

# Get the data cleaned up and properly organized
### Floral survey data
fl_survey <- read_delim("./raw_data/data_plantcover_forDoug_2020.01.17.csv", delim = "\t", 
                           locale = locale(decimal_mark = ",", grouping_mark = ".")) %>%
  na.omit() %>%
  separate(col = species_plant, into = c("genus", "species", "w", "x", "y", "z"), sep = " ") %>%
  unite(species_plant, genus, species, sep = " ") %>%
  select(-c(w,x,y,z)) %>%
  mutate(species_plant = str_replace_all(species_plant, 
                                         c("Arabis caerula" = "Arabis caerulea", # correct nomenclature typos and synonyms
                                           "Arctostaphylos alpina" = "Arctous alpina",
                                           "Calluna vulagris" = "Calluna vulgaris",
                                           "Deschampsia caespitosa" = "Deschampsia cespitosa",
                                           "Gymnadenia bifolia" = "Platanthera bifolia"))) %>%
  separate(species_plant, c("genus", "species"), remove = FALSE) %>%
  dplyr::select(year, dayofyear, site = site_name, snowcover, 
                plant.sp = species_plant, flower.cover = cover_all_flowers) %>%
  mutate(temp.date = as_date(dayofyear),
         day = day(temp.date),
         month = month(temp.date),
         date = make_date(year, month, day))
  
write_csv(fl_survey, "./processed_data/floral_survey.csv")

### Visitation data
network <- read_delim("./raw_data/data_bumblebee_forDoug_2020.04.17.csv", delim = ";", 
                           locale = locale(decimal_mark = ",", grouping_mark = ".")) %>%
  separate(col = plant_sp_latin, into = c("genus", "species", "w", "x", "y", "z"), sep = " ") %>%
  unite(plant_sp_latin, genus, species, sep = " ") %>%
  dplyr::select(-c(w,x,y,z)) %>%
  separate(plant_sp_latin, c("genus", "species"), remove = FALSE) %>%
  dplyr::select(year, dayofyear, site = site_name, trap.time = trapping_time, caste, pollen, 
                bb.sp = bumb_species, bb.sp.lat = bumb_sp_latin, plant.sp.abb = forage_plant, 
                plant.sp = plant_sp_latin) %>%
  mutate(temp.date = as_date(dayofyear),
         day = day(temp.date),
         month = month(temp.date),
         date = make_date(year, month, day)) %>%
  mutate(bb.sp = str_replace_all(bb.sp, c("luco" = "bss","terr" = "bss", "telu" = "bss", "cryp" = "bss", 
                                          "sylv" = "psit", "quad" = "psit", "camp" = "psit", "flav" = "psit", "psyt" = "psit")),
         bb.sp.lat = str_replace_all(bb.sp.lat, c("Bombus terrestris/lucorum" = "Bombus ss",
                                        "Bombus terrestris" = "Bombus ss",
                                        "Bombus lucorum" = "Bombus ss",
                                        "Bombus terrestris/lucorum/lucorum" = "Bombus ss",
                                        "Bombus cryptarum" = "Bombus ss",
                                        "B. barbutellus, B. bohemicus, B. sylvestris" = "Psithyrus",
                                        "Bombus sylvestris" = "Psithyrus",
                                        "Bombus quadricolor" = "Psithyrus",
                                        "Bombus campestris" = "Psithyrus",
                                        "Bombus flavidus" = "Psithyrus")))

write_csv(network, "./processed_data/network.csv")

### Site data
site_data <- read_tsv("./raw_data/SiteLocationFabriceUpdate2020.tsv",
                      locale = locale(decimal_mark = ",", grouping_mark = ".")) %>%
  dplyr::select(site = ID, elev.class = Altitude, management, temp.mean = mean.temp, elev.mean = Mean_altitude, 
                elev.class = Altitude, transect = Transekt, slope.calc = Inklination_degree_calc, slope.est = Inklination_degree_estim,
                elev.min = Min_altitude, elev.max  = Max_altitude, lat = Lat, lon = Long) %>%
  mutate(elev.class = factor(elev.class, ordered = TRUE, levels = c("unten", "mitte", "oben"))) %>%
  mutate(elev.class2 = case_when(
    elev.mean < 1000 ~ "low",
    elev.mean >= 1000 & elev.mean < 1500 ~ "mid",
    elev.mean >= 1500 ~ "high"
  ))

write_csv(site_data, "./processed_data/site_data.csv")

### Floral trait data
floral_k_type <-  read_delim("./BioFlor_traits/BioFlor_Kugler_classification.csv", delim = ";") %>%
  dplyr::select(taxon = 1, k.type = 2) %>% 
  na.omit() %>%
  mutate(taxon = str_replace_all(taxon, c(" x " = " x"))) %>%
  separate(taxon, c("genus", "species"), sep = " ") %>%
  unite(plant.sp, genus, species, sep = " ", remove = FALSE)
  
write_csv(floral_k_type, "./processed_data/floral_k_type.csv")

### Climate data
climate <- read_delim("./raw_data/climate2.txt", delim = "\t") %>%
  rename(site = PlotID) %>%
  mutate(date = dmy(date),
         year = year(date),
         month = month(date),
         site = str_to_lower(site)) %>%
  mutate(gdd = case_when(
    pred_Tmean_day - 10 > 0 ~ pred_Tmean_day - 10,
    pred_Tmean_day - 10 <= 0 ~ 0)
    ) %>%
  group_by(site, year) %>%
  mutate(gdd.cum = round(cumsum(gdd)))

write_csv(climate, "./processed_data/climate.csv")

### TRY data
try <- read_delim("./raw_data/TRY_13_05_2020.txt", delim = "\t")
