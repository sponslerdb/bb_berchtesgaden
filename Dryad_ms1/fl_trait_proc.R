# Packages
library(tidyverse)

# Floral trait data
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
  add_row(plant.sp = "Taraxacum alpinum", plant.genus = "Taraxacum", 
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
          k.type = "NA", k.type.s = "0", k.type.ss = "0") %>%
  add_row(plant.sp = "Tolpis staticifolia", plant.genus = "Tolpis", # Added from BIOFLOR (not sure why it was missing)
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Teucrium montanum", plant.genus = "Teucrium", # Added from BIOFLOR (not sure why it was missing)
          k.type = "5.1", k.type.s = "5.1", k.type.ss = "5") %>%
  add_row(plant.sp = "Silene viscosa", plant.genus = "Silene", # Added from BIOFLOR (not sure why it was missing)
          k.type = "4.2", k.type.s = "4.2", k.type.ss = "4") %>%
  add_row(plant.sp = "Silene pusilla", plant.genus = "Silene", # Added from BIOFLOR (not sure why it was missing)
          k.type = "4.2", k.type.s = "4.2", k.type.ss = "4") %>%
  add_row(plant.sp = "Nigritella nigra", plant.genus = "Nigritella", # Added from BIOFLOR (not sure why it was missing)
          k.type = "7.1", k.type.s = "7.1", k.type.ss = "7") %>%
  add_row(plant.sp = "Melampyrum pratense", plant.genus = "Melampyrum", # Added from BIOFLOR (not sure why it was missing)
          k.type = "5.3", k.type.s = "5.3", k.type.ss = "5") %>%
  add_row(plant.sp = "Leontopodium alpinum", plant.genus = "Leontopodium", # Added from BIOFLOR (not sure why it was missing)
          k.type = "7.2a", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Leontodon sp.", plant.genus = "Leontodon", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Laserpitium siler", plant.genus = "Laserpitium", # Added from BIOFLOR (not sure why it was missing)
          k.type = "1.2a", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Hieracium sp. 1", plant.genus = "Hieracium", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Hieracium sp. 2", plant.genus = "Hieracium", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Hieracium sp. 3", plant.genus = "Hieracium", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Hieracium sp. 4", plant.genus = "Hieracium", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Hieracium glaucum", plant.genus = "Hieracium", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Hieracium glanduliferum", plant.genus = "Hieracium", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Galium uliginosum", plant.genus = "Galium", # Missing
          k.type = "NA", k.type.s = "NA", k.type.ss = "NA") %>%
  add_row(plant.sp = "Erigeron glabratus", plant.genus = "Erigeron", # Unanimous at genus level
          k.type = "7.2c", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Crepis sp.", plant.genus = "Crepis", # Unanimous at genus level
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Clematis alpina", plant.genus = "Clematis", # Added from BIOFLOR (not sure why it was missing)
          k.type = "1.1", k.type.s = "1.1", k.type.ss = "1") %>%
  add_row(plant.sp = "Cicerbita alpina", plant.genus = "Cicerbita", # Added from BIOFLOR (not sure why it was missing)
          k.type = "7.2b", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Chaerophyllum hirsutum ssp. Villarsii", plant.genus = "Chaerophyllum", # Added from BIOFLOR (not sure why it was missing)
          k.type = "1.2a", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Aquilegia vulgaris", plant.genus = "Aquilegia", # Added from BIOFLOR (not sure why it was missing)
          k.type = "3.2", k.type.s = "3.2", k.type.ss = "3") %>%
  add_row(plant.sp = "Anthriscus nitidus", plant.genus = "Anthriscus", # Synomnym
          k.type = "1.2a", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Antennaria carpatica", plant.genus = "Antennaria", # Added from BIOFLOR (not sure why it was missing)
          k.type = "7.2a", k.type.s = "7.2", k.type.ss = "7") %>%
  add_row(plant.sp = "Allium lusitanicum", plant.genus = "Allium", # Added from BIOFLOR (not sure why it was missing)
          k.type = "1.2bd", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Rheum rhabarbarum", plant.genus = "Rheum", # Added from BIOFLOR (not sure why it was missing)
          k.type = "1.2bb", k.type.s = "1.2", k.type.ss = "1") %>%
  add_row(plant.sp = "Mentha aquatica", plant.genus = "Mentha", # Added from BIOFLOR (not sure why it was missing)
          k.type = "2.2", k.type.s = "2.2", k.type.ss = "2") %>%
  add_row(plant.sp = "Gentianella ciliata", plant.genus = "Gentianella", # Added from BIOFLOR (not sure why it was missing)
          k.type = "2.1", k.type.s = "2.1", k.type.ss = "2") %>%
  add_row(plant.sp = "Phyteuma nigrum", plant.genus = "Phyteuma", # synonym
          k.type = "7.1", k.type.s = "7.1", k.type.ss = "7") %>%
  add_row(plant.sp = "Gentianella aspera", plant.genus = "Gentianella", # synonym
          k.type = "2.1", k.type.s = "2.1", k.type.ss = "2") 