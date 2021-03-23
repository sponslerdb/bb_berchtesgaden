library(tidyverse)
library(lubridate)

# Get the data cleaned up and properly organized
### Floral survey data
fl_survey <- read_delim("./data/raw_data/data_plantcover_forDoug_2020.01.17.csv", delim = "\t", 
                           locale = locale(decimal_mark = ",", grouping_mark = ".")) %>%
  na.omit() %>%
  separate(col = species_plant, into = c("genus", "species", "w", "x", "y", "z"), sep = " ") %>%
  unite(species_plant, genus, species, sep = " ") %>%
  dplyr::select(-c(w,x,y,z)) %>%
  mutate(species_plant = str_replace_all(species_plant, 
                                         c("Arabis caerula" = "Arabis caerulea", # correct nomenclature typos and synonyms
                                           "Arctostaphylos alpina" = "Arctous alpina",
                                           "Calluna vulagris" = "Calluna vulgaris",
                                           "Deschampsia caespitosa" = "Deschampsia cespitosa",
                                           "Gymnadenia bifolia" = "Platanthera bifolia",
                                           "Gentiana aspera" = "Gentianella aspera",
                                           "Gentiana ciliata" = "Gentianella ciliata"))) %>%
  separate(species_plant, c("genus", "species"), remove = FALSE) %>%
  dplyr::select(year, dayofyear, site = site_name, snowcover, 
                plant.sp = species_plant, plant.genus = genus, flower.cover = cover_all_flowers) %>%
  mutate(temp.date = as_date(dayofyear),
         day = day(temp.date),
         month = month(temp.date),
         date = make_date(year, month, day)) %>%
  filter(plant.genus != "Huperzia") %>% # because I don't need the clubmoss
  dplyr::select(-temp.date)
  
write_csv(fl_survey, "./data/processed_data/floral_survey.csv")

# 
gentian_problem <- fl_survey %>%
  select(plant.genus, plant.sp) %>%
  filter(plant.genus %in% c("Gentiana", "Gentianella", "Gentianopsis")) %>%
  distinct()

### Visitation data
network <- read_delim("./data/raw_data/data_bumblebee_forDoug_2020.04.17.csv", delim = ";", 
                           locale = locale(decimal_mark = ",", grouping_mark = ".")) %>%
  mutate(plant_sp_latin = str_replace_all(plant_sp_latin, 
                                         c("Arabis caerula" = "Arabis caerulea", # correct nomenclature typos and synonyms
                                           "Arctostaphylos alpina" = "Arctous alpina",
                                           "Calluna vulagris" = "Calluna vulgaris",
                                           "Deschampsia caespitosa" = "Deschampsia cespitosa",
                                           "Gymnadenia bifolia" = "Platanthera bifolia",
                                           "Gentiana aspera" = "Gentianella aspera",
                                           "Gentiana ciliata" = "Gentianella ciliata")))%>%
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
  filter(bb.sp != "cryp") %>% # this species was observed only oncne in three years
  mutate(bb.sp = str_replace_all(bb.sp, c("luco" = "telu", "terr" = "telu", "sylv" = "psit", 
                                          "quad" = "psit", "camp" = "psit", "flav" = "psit", 
                                          "psyt" = "psit")),
         bb.sp.lat = str_replace_all(bb.sp.lat, c("Bombus terrestris/lucorum" = "Bombus ss",
                                        "Bombus terrestris" = "Bombus telu",
                                        "Bombus lucorum" = "Bombus telu",
                                        "Bombus terrestris/lucorum/lucorum" = "Bombus telu",
                                        "B. barbutellus, B. bohemicus, B. sylvestris" = "Psithyrus",
                                        "Bombus sylvestris" = "Psithyrus",
                                        "Bombus quadricolor" = "Psithyrus",
                                        "Bombus campestris" = "Psithyrus",
                                        "Bombus flavidus" = "Psithyrus"))) %>%
  dplyr::select(-temp.date)

write_csv(network, "./data/processed_data/network.csv")

gentian_problem_visitation <- network %>%
  select(plant.genus, plant.sp) %>%
  filter(plant.genus %in% c("Gentiana", "Gentianella", "Gentianopsis")) %>%
  distinct()

### Floral taxonomy
Sys.setenv(ENTREZ_KEY = "ab9a55ec842df6f86a750929aefc69143608")

fl_sp <- fl_survey %>%
  select(plant.sp) %>%
  unique() %>%
  bind_rows(select(network, plant.sp)) %>%
  unique()

fl_tax <- taxize::tax_name(fl_sp$plant.sp, get = "family", db = "ncbi") %>%
  select(plant.sp = query, plant.family = family) %>%
  separate(plant.sp, c("plant.genus", "sp"), remove = FALSE) %>%
  select(plant.sp, plant.genus, plant.family) %>%
  unique()

fl_tax_gapfill <- fl_tax %>% # manually fill in some gaps in the NCBI database
  mutate(plant.family = if_else(is.na(plant.family), 
                                case_when(
                                  plant.genus == "Agrostis" ~ "Poaceae",
                                  plant.genus == "Calamintha" ~ "Lamiaceae",
                                  plant.genus == "Crepis" ~ "Asteraceae",
                                  plant.genus == "Crocus" ~ "Iridaceae",
                                  plant.genus == "Mentha" ~ "Lamiaceae",
                                  plant.genus == "Minuartia" ~ "Caryophyllaceae",
                                  plant.genus == "Papaver" ~ "Papaveraceae",
                                  plant.genus == "Rheum" ~ "Polygonaceae",
                                  plant.genus == "Salix" ~ "Salicaceae",
                                  plant.genus == "Senecio" ~ "Asteraceae"),
                                plant.family)
  )
  
write_csv(fl_tax_gapfill, "./processed_data/floral_tax.csv")
fl_tax_gapfill <- read.csv("./processed_data/floral_tax.csv")

### Site data
site_data <- read_tsv("./raw_data/SiteLocationFabriceUpdate2020.tsv",
                      locale = locale(decimal_mark = ",", grouping_mark = ".")) %>%
  dplyr::select(site = ID, elev.class = Altitude, management, temp.mean = mean.temp, elev.mean = Mean_altitude, 
                elev.class = Altitude, transect = Transekt, slope.calc = Inklination_degree_calc, slope.est = Inklination_degree_estim,
                elev.min = Min_altitude, elev.max  = Max_altitude, lat = Lat, lon = Long) %>%
  mutate(elev.class = factor(elev.class, ordered = TRUE, levels = c("unten", "mitte", "oben"))) %>%
  mutate(elev.class2 = case_when(
    elev.mean <= 850 ~ "low",
    elev.mean > 850 & elev.mean <= 1500 ~ "mid",
    elev.mean > 1500 ~ "high"
  ))

write_csv(site_data, "./processed_data/site_data.csv")

### Floral trait data
floral_k_type <-  read_delim("./data/raw_data/BioFlor_Kugler_classification.csv", delim = ";") %>%
  dplyr::select(taxon = 1, k.type = 2) %>% 
  na.omit() %>%
  mutate(taxon = str_replace_all(taxon, c(" x " = " x"))) %>%
  separate(taxon, c("genus", "species"), sep = " ") %>%
  unite(plant.sp, genus, species, sep = " ", remove = FALSE) %>%
  select(plant.sp, plant.genus = genus, k.type) %>%
  mutate(k.type.s = str_extract(k.type, "[-+]?[0-9]*\\.?[0-9]+"), # make simplified (and more simplified) k types
         k.type.ss = str_extract(k.type, "[-+]?[0-9]*")) %>%
  
  # Fill in ktype data for plants present in my samples but missing from the Bioflor database
  add_row(plant.sp = "Phyteuma orbiculare", plant.genus = "Phyteuma", 
          k.type = "NA", k.type.s = "7.1", k.type.ss = "7") %>% # Following Neumayer and Paulus (1999)
  add_row(plant.sp = "Phyteuma spicatum", plant.genus = "Phyteuma",
          k.type = "NA", k.type.s = "7.1", k.type.ss = "7") %>% # Following Neumayer and Paulus (1999)
  add_row(plant.sp = "Aconitum vulparia", plant.genus = "Aconitum",
          k.type = "NA", k.type.s = "5.2", k.type.ss = "5") %>% # All scored Aconitum were classified as 5.2
  add_row(plant.sp = "Cardamine enneaphyllos", plant.genus = "Cardamine", 
          k.type = "NA", k.type.s = "1.2", k.type.ss = "1") %>% # All scored Cardamine were classified as 1.2
  add_row(plant.sp = "Euphrasia picta", plant.genus = "Euphrasia", 
          k.type = "NA", k.type.s = "5.1", k.type.ss = "5") %>% # All scored Euphrasia were classified as 5.1
  add_row(plant.sp = "Euphrasia rostkoviana", plant.genus = "Euphrasia", 
          k.type = "NA", k.type.s = "5.1", k.type.ss = "5") %>% # All scored Euphrasia were classified as 5.1
  add_row(plant.sp = "Lamium galeobdolon", plant.genus = "Lamium", 
          k.type = "NA", k.type.s = "5.1", k.type.ss = "5") %>% # All scored Lamium were classified as 5.1
  add_row(plant.sp = "Mentha aquatilis", plant.genus = "Mentha", 
          k.type = "NA", k.type.s = "2.2", k.type.ss = "2") %>% # All scored Mentha were classified as 2.2
  add_row(plant.sp = "Rumex alpestris", plant.genus = "Rumex", 
          k.type = "NA", k.type.s = "0.0", k.type.ss = "0") %>% # All scored Rumex were classified as 0
  add_row(plant.sp = "Rheum barbarum", plant.genus = "Rheum", 
          k.type = "NA", k.type.s = "7.1", k.type.ss = "7") %>% # My fiat declaration
  add_row(plant.sp = "Salix sp", plant.genus = "Salix", 
          k.type = "NA", k.type.s = "9.0", k.type.ss = "9") %>% # All scored Salix were classified as 9
  add_row(plant.sp = "Salix spec.", plant.genus = "Salix", 
          k.type = "NA", k.type.s = "9.0", k.type.ss = "9") %>% # All scored Salix were classified as 9
  add_row(plant.sp = "Silene flos-cuculi", plant.genus = "Silene", 
          k.type = "NA", k.type.s = "4.2", k.type.ss = "4") %>% # All scored Silene were classified as 4.2
  add_row(plant.sp = "Stachys alopecuros", plant.genus = "Stachys", 
          k.type = "NA", k.type.s = "5.1", k.type.ss = "5") %>% # All scored Stachys were classified as 5.1
  add_row(plant.sp = "Stachys officinalis", plant.genus = "Stachys", 
          k.type = "NA", k.type.s = "5.1", k.type.ss = "5") %>% # All scored Stachys were classified as 5.1
  add_row(plant.sp = "Taraxacum officinale", plant.genus = "Taraxacum", 
          k.type = "NA", k.type.s = "7.2", k.type.ss = "7") %>% # All scored Taraxacum (treated as a species group by Kugler) belong to 7.2
  add_row(plant.sp = "Alchemilla conjuncta", plant.genus = "Alchemilla", # All scored Alchemilla were classified as 1.2 
          k.type = "NA", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Arabis pumila", plant.genus = "Arabis", # All scored Arabis were classified as 1.2 
          k.type = "NA", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Arctous alpina", plant.genus = "Arctous", # All scored Ericaceae were classified as 3.1, and A. alpina is clearly a bell flower 
          k.type = "NA", k.type.s = "3.1", k.type.ss = "3") %>%
  add_row(plant.sp = "Calamintha clinopodium", plant.genus = "Calamintha", # Taxonomic synonym Clinopodium nepeta 
          k.type = "NA", k.type.s = "5.1", k.type.ss = "5") %>%
  add_row(plant.sp = "Carex leporina", plant.genus = "Carex", # All scored Carex were classified as 0 
          k.type = "NA", k.type.s = "0.0", k.type.ss = "0") %>%
  add_row(plant.sp = "Ficaria verna", plant.genus = "Ficaria", # Taxonomic synonym Ranunculus ficaria
          k.type = "NA", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Huperzia selago", plant.genus = "Huperzia", # It's a fern
          k.type = "NA", k.type.s = "NA", k.type.ss = "NA") %>%
  add_row(plant.sp = "Minuartia gerardii", plant.genus = "Minuartia", # All scored Minuartia were classified as 1.2 
          k.type = "NA", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Papaver aurantiacum", plant.genus = "Papaver", # All Papaver are 1.1 
          k.type = "NA", k.type.s = "1.1", k.type.ss = "1") %>%
  add_row(plant.sp = "Polygonum viviparum", plant.genus = "Polygonum", # Synonym: Bistorta vivipara
          k.type = "NA", k.type.s = "3.2", k.type.ss = "3") %>%
  add_row(plant.sp = "Sesleria caerulea", plant.genus = "Sesleria", # It's a grass
          k.type = "NA", k.type.s = "0", k.type.ss = "0")
  

write_csv(floral_k_type, "./data/processed_data/floral_k_type.csv")

gentian_problem_traits <- floral_k_type %>%
  select(plant.genus, plant.sp, k.type.ss) %>%
  filter(plant.genus %in% c("Gentiana", "Gentianella", "Gentianopsis")) %>%
  distinct()

floral_m_type <-  read_delim("./BioFlor_traits/BioFlor_Mueller_classification.csv", delim = ";") %>%
  dplyr::select(taxon = 1, m.type = 2) %>% 
  na.omit() %>%
  mutate(taxon = str_replace_all(taxon, c(" x " = " x"))) %>%
  separate(taxon, c("genus", "species"), sep = " ") %>%
  unite(plant.sp, genus, species, sep = " ", remove = FALSE) %>%
  select(plant.sp, plant.genus = genus, m.type) 

floral_color <-  read_delim("./BioFlor_traits/BioFlor_flower_color.csv", delim = ";") %>%
  slice(-1) %>%
  select(taxon = 1, color = 2) %>% 
  na.omit() %>%
  head(-1) %>%
  mutate(taxon = str_replace_all(taxon, c(" x " = " x"))) %>%
  separate(taxon, c("genus", "species"), sep = " ") %>%
  unite(plant.sp, genus, species, sep = " ", remove = FALSE) %>%
  select(plant.sp, plant.genus = genus, color) %>%
  
  # Fill in color data for plants present in my samples but missing from the Bioflor database
  # Colors based on review of iNaturalist "research-grade" observations or Wikimedia commons. Where color is too variable to classify, "various_colors". 
  # Where essentially absent as in conifers, "NA".
  add_row(plant.sp = "Euphrasia rostkoviana", plant.genus = "Euphrasia", color = "var") %>%
  add_row(plant.sp = "Euphrasia officinalis", plant.genus = "Euphrasia", color = "var") %>%
  add_row(plant.sp = "Taraxacum officinale", plant.genus = "Taraxacum", color = "yellow") %>%
  add_row(plant.sp = "Juniperus communis", plant.genus = "Juniperus", color = "green") %>%
  add_row(plant.sp = "Gentiana aspera", plant.genus = "Gentiana", color = "purple") %>%
  add_row(plant.sp = "Gentianella aspera", plant.genus = "Gentianella", color = "purple") %>%
  add_row(plant.sp = "Carex ericetorum", plant.genus = "Carex", color = "yellow") %>%
  add_row(plant.sp = "Euphrasia picta", plant.genus = "Euphrasia", color = "various colors") %>%
  add_row(plant.sp = "Gentiana ciliata", plant.genus = "Gentiana", color = "violet") %>%
  add_row(plant.sp = "Gentianopsis ciliata", plant.genus = "Gentiana", color = "violet") %>%
  add_row(plant.sp = "Stachys alopecuros", plant.genus = "Stachys", color = "yellow") %>%
  add_row(plant.sp = "Betonica alopecuros", plant.genus = "Betonica", color = "yellow") %>%
  add_row(plant.sp = "Cardamine enneaphyllos", plant.genus = "Cardamine", color = "yellow") %>%
  add_row(plant.sp = "Lamium galeobdolon", plant.genus = "Lamium", color = "yellow") %>%
  add_row(plant.sp = "Larix decidua", plant.genus = "Larix", color = "red") %>%
  add_row(plant.sp = "Silene flos-cuculi", plant.genus = "Silene", color = "pink") %>%
  add_row(plant.sp = "Aconitum vulparia", plant.genus = "Aconitum", color = "yellow") %>%
  add_row(plant.sp = "Rumex alpestris", plant.genus = "Rumex", color = "red") %>%
  add_row(plant.sp = "Mentha aquatilis", plant.genus = "Mentha", color = "lilac") %>%
  add_row(plant.sp = "Mentha aquatica", plant.genus = "Mentha", color = "lilac") %>%
  add_row(plant.sp = "Salix spec.", plant.genus = "Salix", color = "var") %>%
  add_row(plant.sp = "Stachys officinalis", plant.genus = "Stachys", color = "pink") %>%
  add_row(plant.sp = "Betonica officinalis", plant.genus = "Betonica", color = "pink") %>%
  add_row(plant.sp = "Rheum barbarum", plant.genus = "Rheum", color = "yellow") %>%
  add_row(plant.sp = "Rheum rhabarbarum", plant.genus = "Rheum", color = "yellow")
  

write_csv(floral_color, "./processed_data/floral_color.csv")

### Bee colors from FReD database, downloaded December 8, 2020. http://www.reflectance.co.uk//index.php
fred <- read_csv("./raw_data/FreD.csv") %>%
  unite(plant.sp, c(Genus, Species), sep = " ") %>%
  
  # Fill in missing data where possible
  
  add_row(plant.sp = "Erica carnea", BeeColour = "blue", HumanColour = "pink") %>% # taxonomic synonym
  add_row(plant.sp = "Centaurea jacea", BeeColour = "blue-green", HumanColour = "pink") %>% # Raine et al. (2007) PLoS ONE 
  add_row(plant.sp = "Salvia verticillata", BeeColour = "uv-blue", HumanColour = "lilac") %>% # Raine et al. (2007) PLoS ONE
  add_row(plant.sp = "Epilobium angustifolium", BeeColour = "blue", HumanColour = "pink") %>% # Raine et al. (2007) PLoS ONE
  add_row(plant.sp = "Lamium galeobdolon", BeeColour = "uv-green", HumanColour = "yellow") # Raine et al. (2007) PLoS ONE

survey_species <- read_csv("./processed_data/floral_survey.csv") %>%
  select(plant.genus, plant.sp) %>%
  distinct()

visits_species <- read_csv("./processed_data/network.csv") %>%
  select(plant.genus, plant.sp) %>%
  distinct()
  
fred_survey <- left_join(survey_species, fred) %>%
  select(plant.sp, BeeColour) %>%
  distinct()

fred_visits <- left_join(visits_species, fred) %>%
  select(plant.sp, BeeColour) %>%
  distinct()

missing_freds <- bind_rows(fred_survey, fred_visits) %>%
  distinct() %>%
  filter(is.na(BeeColour))

missing_freds.survey <- fred_survey %>%
  distinct() %>%
  filter(is.na(BeeColour))

missing_freds.visits <- fred_visits %>%
  distinct() %>%
  filter(is.na(BeeColour))


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

### Weather data
weather2010 <- read_delim("./raw_data/field-data_bumblebees_2010_raw-data.csv", delim = ";") %>%
  dplyr::select(site = site_name, date, weather, temperature) %>% 
  unique() %>%
  mutate(date = dmy(date))

weather2011 <- read_delim("./raw_data/field-data_bumblebees_2011_raw-data.csv", delim = ";") %>%
  dplyr::select(site = site_name, date, weather, temperature) %>% 
  unique() %>%
  mutate(date = dmy(date))

weather2012 <- read_delim("./raw_data/field-data_bumblebees_2012_raw-data.csv", delim = ";") %>%
  dplyr::select(site = site_name, date, weather, temperature) %>%
  mutate(site = str_replace_all(site, "salix", "")) %>%
  mutate(site = str_replace_all(site, "garden", "")) %>%
  unique() 

weather <- bind_rows(weather2010, weather2011, weather2012) %>%
  write_csv("./processed_data/weather.csv")

temperature <- read_delim("./raw_data/climate2.txt", delim = "\t") %>%
  rename(site = PlotID) %>%
  mutate(date = dmy(date),
         site = str_to_lower(site)) %>%
  write_csv("./processed_data/temperature.csv")
