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
                plant.sp = species_plant, plant.genus = genus, flower.cover = cover_all_flowers) %>%
  mutate(temp.date = as_date(dayofyear),
         day = day(temp.date),
         month = month(temp.date),
         date = make_date(year, month, day))
  
write_csv(fl_survey, "./processed_data/floral_survey.csv")

### Visitation data
network <- read_delim("./raw_data/data_bumblebee_forDoug_2020.04.17.csv", delim = ";", 
                           locale = locale(decimal_mark = ",", grouping_mark = ".")) %>%
  mutate(plant_sp_latin = str_replace_all(plant_sp_latin, 
                                         c("Arabis caerula" = "Arabis caerulea", # correct nomenclature typos and synonyms
                                           "Arctostaphylos alpina" = "Arctous alpina",
                                           "Calluna vulagris" = "Calluna vulgaris",
                                           "Deschampsia caespitosa" = "Deschampsia cespitosa",
                                           "Gymnadenia bifolia" = "Platanthera bifolia")))%>%
  separate(col = plant_sp_latin, into = c("genus", "species", "w", "x", "y", "z"), sep = " ") %>%
  unite(plant_sp_latin, genus, species, sep = " ") %>%
  dplyr::select(-c(w,x,y,z)) %>%
  separate(plant_sp_latin, c("genus", "species"), remove = FALSE) %>%
  dplyr::select(year, dayofyear, site = site_name, trap.time = trapping_time, caste, pollen, 
                bb.sp = bumb_species, bb.sp.lat = bumb_sp_latin, plant.sp.abb = forage_plant, 
                plant.sp = plant_sp_latin, plant.genus = genus) %>%
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

### Floral taxonomy
Sys.setenv(ENTREZ_KEY = "ab9a55ec842df6f86a750929aefc69143608")

fl_sp <- survey %>%
  select(plant.sp) %>%
  unique() %>%
  bind_rows(select(network, plant.sp)) %>%
  unique()

fl_tax <- tax_name(fl_sp$plant.sp, get = "family", db = "ncbi") %>%
  select(plant.sp = query, plant.family = family) %>%
  separate(plant.sp, c("plant.genus", "plant.sp"), remove = FALSE) %>%
  select(plant.genus, plant.family) %>%
  unique()

fl_tax_gapfill <- fl_tax %>% # manually fill in some gaps in the NCBI database
  mutate(plant.family = case_when(
    plant.genus == "Crocus" ~ "Iridaceae",
    plant.genus == "Rheum" ~ "Polygonaceae",
    plant.genus == "Mentha" ~ "Lamiaceae",
    plant.genus %in% c("Senecio", "Crepis") ~ "Asteraceae",
    plant.genus == "Salix" ~ "Salicaceae",
    !plant.genus %in% c("Senecio", "Crocus", "Crepis", "Mentha", "Salix", "Rheum") ~ plant.family
  ))

write_csv(fl_tax_gapfill, "./processed_data/floral_tax.csv")


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
    elev.mean >= 1750 ~ "high"
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

### Parasite data
#### 2010
parasites2010 <- read_delim("./raw_data/lab-data_parasites_2010_raw-data.csv", delim = ";") %>%
  select(-c(year, ovary_development, fatbody_look, comment, dissector)) %>%
  gather(var, value, -sample_number) %>%
  mutate(value = as.numeric(value),
         var = str_replace_all(var, "-", "_")) %>%
  na.omit() %>%
  spread(var, value) %>%
  filter(!is.na(size)) %>% #remove the few cases where size was not measured 
  replace(is.na(.), 0) %>%
  mutate(ecto.mites = ecto_mites_l + ecto_mites_m + ecto_mites_s + ecto_mites_xs,
         conopids = parasitoid_conopidae_dead_larvae +
           parasitoid_conopidae_dead_larvae +
           parasitoid_conopidae_dorsal_egg +
           parasitoid_conopidae_dorsal_larvae +
           parasitoid_conopidae_ventral_egg +
           parasitoid_conopidae_ventral_larvae) %>%
  select(sample_number, size, ecto.mites, conopids, 
         nosema = parasite_nosema, apicystis = parasite_apicystis,
         crithidia = parasite_crithidia, braconids = parasitoid_braconidae_larvae,
         tracheal.mites = tracheal_mites, nematodes = nematode) %>%
  gather(var, value, -sample_number, -size) %>%
  mutate(value = as.numeric(value),
         bin.value = case_when(
           value > 0 ~ TRUE,
           value == 0 ~ FALSE
         )) %>%
  na.omit() 
  
sampling2010 <- read_delim("./raw_data/field-data_bumblebees_2010_raw-data.csv", delim = ";") %>%
  select(site = site_name, elev.mean = site_elevation, year, date, yday = dayofyear, bb.sp = bumb_species, 
         caste, plant.sp = forage_plant, sample_number) %>%
  mutate(date = dmy(date),
         caste = str_trim(caste))

psite_2010 <- parasites2010 %>%
  left_join(sampling2010, by = c("sample_number")) %>%
  select(site, date, year, yday, elev.mean, bb.sp, size, caste, 
         plant.sp, var, value, bin.value, sample_number)

#### 2011
parasites2011 <- read_delim("./raw_data/lab-data_parasites_2011_raw-data.csv", delim = ";") %>%
  select(-c(year, ovary_development, fatbody_look, comment, dissector)) %>%
  gather(var, value, -sample_number) %>%
  mutate(value = as.numeric(value),
         var = str_replace_all(var, "-", "_")) %>%
  na.omit() %>%
  spread(var, value) %>%
  filter(!is.na(size)) %>% #remove the few cases where size was not measured 
  replace(is.na(.), 0) %>%
  mutate(ecto.mites = ecto_mites_l + ecto_mites_m + ecto_mites_s + ecto_mites_xs,
         conopids = parasitoid_conopidae_dead_larvae +
           parasitoid_conopidae_dead_larvae +
           parasitoid_conopidae_dorsal_egg +
           parasitoid_conopidae_dorsal_larvae +
           parasitoid_conopidae_ventral_egg +
           parasitoid_conopidae_ventral_larvae) %>%
  select(sample_number, size, ecto.mites, conopids, 
         braconids = parasitoid_braconidae_larvae,
         tracheal.mites = tracheal_mites, nematodes = nematode) %>%
  gather(var, value, -sample_number, -size) %>%
  mutate(value = as.numeric(value),
         bin.value = case_when(
           value > 0 ~ TRUE,
           value == 0 ~ FALSE
         )) %>%
  na.omit()  %>%
  mutate(sample_number = as.character(sample_number))

sampling2011 <- read_delim("./raw_data/field-data_bumblebees_2011_raw-data.csv", delim = ";") %>%
  select(site = site_name, elev.mean = site_elevation, year, date, yday = dayofyear, bb.sp = bumb_species, 
         caste, plant.sp = forage_plant, sample_number) %>%
  mutate(date = dmy(date),
         caste = str_trim(caste)) %>%
  mutate(sample_number = as.character(sample_number))

psite_2011 <- parasites2011 %>%
  left_join(sampling2011, by = c("sample_number")) %>%
  select(site, date, year, yday, elev.mean, bb.sp, size, caste, 
         plant.sp, var, value, bin.value, sample_number)

ggplot(psite_2011, aes(bin.value)) +
  geom_bar() +
  facet_wrap(~var, scale = "free")

ggplot(psite_2011, aes(value)) +
  geom_bar() +
  facet_wrap(~var, scale = "free")

#### 2012
parasites2012 <- read_delim("./raw_data/lab-data_parasites_2012_raw-data_mod.csv", delim = "\t") %>%
  select(-c(year, ovary_development, comment, dissector)) %>%
  gather(var, value, -sample_number) %>%
  mutate(value = as.numeric(value),
         var = str_replace_all(var, "-", "_")) %>%
  na.omit() %>%
  spread(var, value) %>%
  filter(!is.na(size)) %>% #remove the few cases where size was not measured 
  replace(is.na(.), 0) %>%
  mutate(ecto.mites = ecto_mites_l + ecto_mites_m + ecto_mites_s + ecto_mites_xs,
         conopids = parasitoid_conopidae_dead_larvae +
           parasitoid_conopidae_dead_larvae +
           parasitoid_conopidae_dorsal_egg +
           parasitoid_conopidae_dorsal_larvae +
           parasitoid_conopidae_ventral_egg +
           parasitoid_conopidae_ventral_larvae) %>%
  select(sample_number, size, ecto.mites, conopids, 
         nosema = parasite_nosema, apicystis = parasite_apicystis,
         crithidia = parasite_crithidia, braconids = parasitoid_braconidae_larvae,
         tracheal.mites = tracheal_mites, nematodes = nematode) %>%
  gather(var, value, -sample_number, -size) %>%
  mutate(value = as.numeric(value),
         bin.value = case_when(
           value > 0 ~ TRUE,
           value == 0 ~ FALSE
         )) %>%
  na.omit() %>%
  mutate(sample_number = as.character(sample_number))

sampling2012 <- read_delim("./raw_data/field-data_bumblebees_2012_raw-data.csv", delim = ";") %>%
  select(site = site_name, elev.mean = site_elevation, year, date, yday = dayofyear, bb.sp = bumb_species, 
         caste, plant.sp = forage_plant, sample_number) %>%
  mutate(caste = str_trim(caste)) %>%
  filter(sample_number > 0) %>%
  mutate(sample_number = as.character(sample_number))

psite_2012 <- parasites2012 %>%
  left_join(sampling2012, by = c("sample_number")) %>%
  select(site, date, year, yday, elev.mean, bb.sp, size, caste, 
         plant.sp, var, value, bin.value, sample_number)

parasites <- bind_rows(psite_2010, psite_2011, psite_2012)

write_csv(parasites, "./processed_data/parasites.csv")
