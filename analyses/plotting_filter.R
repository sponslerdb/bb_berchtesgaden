### Need to run outside main script to avoid memory crashes

library(tidyverse)
library(lubridate)

site_data <- read_csv("./processed_data/site_data.csv") %>%
  mutate(elev.class = factor(elev.class, levels = c("oben", "mitte", "unten")), # rename in English
         elev.class2 = factor(elev.class2, levels = c("high", "mid", "low")))

net <- read_csv("./processed_data/network.csv") %>%
  dplyr::select(site, date, bb.sp, plant.sp) # drop unused columns

survey <- read_csv("./processed_data/floral_survey.csv") %>%
  semi_join(net, by = "plant.sp") %>% # Consider only visited plant species (anytime anywhere, not specifically for a given site*date)
  dplyr::select(site, date, plant.sp, flower.cover) %>%
  group_by(plant.sp) %>%
  mutate(present = case_when(
    max(flower.cover) > 0 ~ TRUE,
    max(flower.cover) == 0 ~ FALSE
  )) %>%
  filter(present == TRUE) 

### This could be so much more beautiful, but I don't have it in me right now. This works.
seq_yday <- function(x){
  seq((x - 10), (x + 10), 0.5)
}

seq_elev <- function(x){
  seq((x - 150), (x + 150), 1)
}

df <- function(x) {
  as.data.frame(x)
}

survey_range.2010 <- survey %>%
  ungroup() %>%
  left_join(site_data) %>%
  mutate(yday = yday(date),
         year = year(date)) %>%
  filter(year == 2010) %>%
  select(yday, elev.mean) %>%
  unique() %>%
  arrange(yday) %>%
  mutate(ydays = map(yday, seq_yday),
         elevs = map(elev.mean, seq_elev)) %>%
  mutate(ydays = map(ydays, df),
         elevs = map(elevs, df)) %>%
  unnest(cols = ydays) %>%
  select(yday = x, elevs) %>%
  unnest(cols = elevs) %>%
  unique() %>%
  select(yday, elev.mean = x) %>%
  mutate(year = rep(2010, n()))

survey_range.2011 <- survey %>%
  ungroup() %>%
  left_join(site_data) %>%
  mutate(yday = yday(date),
         year = year(date)) %>%
  filter(year == 2011) %>%
  select(yday, elev.mean) %>%
  unique() %>%
  arrange(yday) %>%
  mutate(ydays = map(yday, seq_yday),
         elevs = map(elev.mean, seq_elev)) %>%
  mutate(ydays = map(ydays, df),
         elevs = map(elevs, df)) %>%
  unnest(cols = ydays) %>%
  select(yday = x, elevs) %>%
  unnest(cols = elevs) %>%
  unique() %>%
  select(yday, elev.mean = x) %>%
  mutate(year = rep(2011, n()))

survey_range.2012 <- survey %>%
  ungroup() %>%
  left_join(site_data) %>%
  mutate(yday = yday(date),
         year = year(date)) %>%
  filter(year == 2012) %>%
  select(yday, elev.mean) %>%
  unique() %>%
  arrange(yday) %>%
  mutate(ydays = map(yday, seq_yday),
         elevs = map(elev.mean, seq_elev)) %>%
  mutate(ydays = map(ydays, df),
         elevs = map(elevs, df)) %>%
  unnest(cols = ydays) %>%
  select(yday = x, elevs) %>%
  unnest(cols = elevs) %>%
  unique() %>%
  select(yday, elev.mean = x) %>%
  mutate(year = rep(2012, n()))

bb_range.2010 <- net %>%
  ungroup() %>%
  left_join(site_data) %>%
  mutate(yday = yday(date),
         year = year(date)) %>%
  filter(year == 2010) %>%
  select(yday, elev.mean) %>%
  unique() %>%
  arrange(yday) %>%
  mutate(ydays = map(yday, seq_yday),
         elevs = map(elev.mean, seq_elev)) %>%
  mutate(ydays = map(ydays, df),
         elevs = map(elevs, df)) %>%
  unnest(cols = ydays) %>%
  select(yday = x, elevs) %>%
  unnest(cols = elevs) %>%
  unique() %>%
  select(yday, elev.mean = x) %>%
  mutate(year = rep(2010, n()))

bb_range.2011 <- net %>%
  ungroup() %>%
  left_join(site_data) %>%
  mutate(yday = yday(date),
         year = year(date)) %>%
  filter(year == 2011) %>%
  select(yday, elev.mean) %>%
  unique() %>%
  arrange(yday) %>%
  mutate(ydays = map(yday, seq_yday),
         elevs = map(elev.mean, seq_elev)) %>%
  mutate(ydays = map(ydays, df),
         elevs = map(elevs, df)) %>%
  unnest(cols = ydays) %>%
  select(yday = x, elevs) %>%
  unnest(cols = elevs) %>%
  unique() %>%
  select(yday, elev.mean = x) %>%
  mutate(year = rep(2011, n()))

bb_range.2012 <- net %>%
  ungroup() %>%
  left_join(site_data) %>%
  mutate(yday = yday(date),
         year = year(date)) %>%
  filter(year == 2012) %>%
  select(yday, elev.mean) %>%
  unique() %>%
  arrange(yday) %>%
  mutate(ydays = map(yday, seq_yday),
         elevs = map(elev.mean, seq_elev)) %>%
  mutate(ydays = map(ydays, df),
         elevs = map(elevs, df)) %>%
  unnest(cols = ydays) %>%
  select(yday = x, elevs) %>%
  unnest(cols = elevs) %>%
  unique() %>%
  select(yday, elev.mean = x) %>%
  mutate(year = rep(2012, n()))

survey_range <- bind_rows(survey_range.2010, survey_range.2011, survey_range.2012)
bb_range <- bind_rows(bb_range.2010, bb_range.2011, bb_range.2012)
