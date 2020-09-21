
boot_net_sp_FL <- function(net, n, iter) {
  
  out <- data.frame()
  
  for(i in 1:iter) {
    temp <- net %>%
      dplyr::select(site, date) %>%
      unique() %>%
      group_by(site) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
      summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
      group_by(site) %>%
      filter(length(unique(bb.sp)) > 1 & length(unique(plant.sp)) > 1) %>%
      ungroup() %>%
      filter(freq > 1) %>% # remove singletons
      dplyr::select(higher = bb.sp, lower = plant.sp, webID = site, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      map(specieslevel, level = "lower", index = c("d", "species strength", "proportional generality")) %>%
      map(rownames_to_column) %>%
      map(as_tibble) %>%
      bind_rows(.id = "site") %>% # Collapse list into big data frame
      rename(sp = rowname) %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out %>%
    gather(metric, value, -c(sp, site, iteration)) %>% # get mean value across iterations for each metric
    group_by(sp, site, metric) %>%
    summarize(value = mean(value)) %>%
    left_join(site_data) %>%
    left_join(fl_traits, by = c("sp" = "plant.sp"))
}


#### Rarefied network analysis by colony phase
```{r, include=FALSE}
# network level version 2
boot_net_nt2 <- function(net, n, iter) {
  
  out <- data.frame()
  
  t_filter <- net %>%
    dplyr::select(site, date, phase, year) %>%
    unique() %>%
    group_by(site, phase, year) %>%
    mutate(transects = n()) %>%
    filter(transects >= n)
  
  for(i in 1:iter) {
    temp <- net %>%
      semi_join(t_filter) %>%
      dplyr::select(site, date, phase, year) %>%
      unique() %>%
      group_by(site, phase, year) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, phase, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
      summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
      #filter(freq > 1) %>%
      unite(webID, c(site, phase, year), sep = "_") %>%
      dplyr::select(higher = bb.sp, lower = plant.sp, webID, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      map(networklevel, weighted = TRUE, level = "higher",
          index = c("web asymmetry", "weighted NODF", "H2", "C score", "niche overlap", "generality")) %>% 
      data.frame() %>%
      t() %>%
      data.frame() %>%
      rownames_to_column() %>%
      dplyr::select(id = rowname, everything()) %>%
      as_tibble() %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out  %>%
    gather(metric, value, -c(id, iteration)) %>% # get mean value across iterations for each metric
    group_by(id, metric) %>%
    summarize(value = mean(value)) %>%
    separate(id, sep = "_", into = c("site", "phase", "year")) %>%
    mutate(phase = factor(phase, levels = c("founding", "buildup", "reproductive"))) %>%
    left_join(site_data)
}
# species level version 2
boot_net_sp2 <- function(net, n, iter) {
  
  out <- data.frame()
  
  t_filter <- net %>%
    dplyr::select(site, date, phase, year) %>%
    unique() %>%
    group_by(site, phase, year) %>%
    mutate(transects = n()) %>%
    filter(transects >= n)
  
  for(i in 1:iter) {
    temp <- net %>%
      semi_join(t_filter) %>%
      dplyr::select(site, date, phase, year) %>%
      unique() %>%
      group_by(site, phase, year) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, phase, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
      summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
      #filter(freq > 1) %>%
      unite(webID, c(site, phase, year), sep = "_") %>%
      dplyr::select(higher = bb.sp, lower = plant.sp, webID, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      map(specieslevel, level = "higher", index = c("d", "species strength")) %>%
      map(rownames_to_column) %>%
      map(as_tibble) %>%
      bind_rows(.id = "id") %>% # Collapse list into big data frame
      rename(sp = rowname) %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out %>%
    gather(metric, value, -c(sp, id, iteration)) %>% # get mean value across iterations for each metric
    group_by(sp, id, metric) %>%
    summarize(value = mean(value)) %>%
    separate(id, sep = "_", into = c("site", "phase", "year")) %>%
    mutate(phase = factor(phase, levels = c("founding", "buildup", "reproductive"))) %>%
    left_join(site_data) %>%
    left_join(bb_traits, by = c("sp" = "bb.sp"))
}
```


### Rarefied network analysis by floral k.type instead of floral species
```{r}
# boot_net function: rarefies networks by iteratively subsampling more frequently sampled networks to match the transect count of the least frequently sampled one
# species level
boot_net_sp_k.type.s <- function(net, n, iter) {
  
  out <- data.frame()
  
  for(i in 1:iter) {
    temp <- net %>%
      dplyr::select(site, date) %>%
      unique() %>%
      group_by(site) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, bb.sp, k.type.s) %>% # group by site, date, bee species, plant species
      summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
      na.omit() %>%
      group_by(site) %>%
      filter(length(unique(bb.sp)) > 1 & length(unique(k.type.s)) > 1) %>%
      ungroup() %>%
      filter(freq > 1) %>% # remove singletons
      dplyr::select(higher = bb.sp, lower = k.type.s, webID = site, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      map(specieslevel, level = "higher", index = c("d", "species strength", "proportional generality")) %>%
      map(rownames_to_column) %>%
      map(as_tibble) %>%
      bind_rows(.id = "site") %>% # Collapse list into big data frame
      rename(sp = rowname) %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out %>%
    gather(metric, value, -c(sp, site, iteration)) %>% # get mean value across iterations for each metric
    group_by(sp, site, metric) %>%
    summarize(value = mean(value)) %>%
    left_join(site_data) %>%
    left_join(bb_traits, by = c("sp" = "bb.sp"))
}

# network level
boot_net_nt_k.type.s <- function(net, n, iter) {
  
  out <- data.frame()
  
  for(i in 1:iter) {
    temp <- net %>%
      dplyr::select(site, date) %>%
      unique() %>%
      group_by(site) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, bb.sp, k.type.s) %>% # group by site, date, bee species, plant species
      summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
      na.omit() %>%
      group_by(site) %>%
      filter(length(unique(bb.sp)) > 1 & length(unique(k.type.s)) > 1) %>%
      ungroup() %>%
      dplyr::select(higher = bb.sp, lower = k.type.s, webID = site, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      map(networklevel, weighted = TRUE, level = "higher", index = c("H2", "weighted NODF", "niche overlap", "C score", "generality", "web asymmetry", "weighted connectance")) %>% # We probably care mainly about the BB level of the network
      data.frame() %>%
      t() %>%
      data.frame() %>%
      rownames_to_column() %>%
      dplyr::select(site = rowname, everything()) %>%
      as_tibble() %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out  %>%
    gather(metric, value, -c(site, iteration)) %>% # get mean value across iterations for each metric
    group_by(site, metric) %>%
    summarize(value = mean(value)) %>%
    left_join(site_data)
}

bb_traits <- read_csv("../processed_data/bb_traits.csv") %>%
  dplyr::select(subgenus, bb.sp, pbl.w, pbl.w.class)

fl_traits <- read.csv("../processed_data/floral_k_type_nom_fix.csv") %>%
  dplyr::select(plant.sp, k.type) %>%
  mutate(k.type.s = str_extract(k.type, "[-+]?[0-9]*\\.?[0-9]+"),
         k.type.ss = str_extract(k.type, "[-+]?[0-9]*"))

fl_tax <- read_csv("../processed_data/floral_tax.csv") 

breaks <- read_csv("./output/breaks.csv") %>%
  mutate(year = factor(year))

```


### Rarefied network analysis by k.type instead of by plant sp
```{r, message = FALSE}
### Create yearly subsets to feed to boot_net
net_2012 <- filter(net, year == "2012")
net_2011 <- filter(net, year == "2011")
net_2010 <- filter(net, year == "2010")

### Species level
spmet_2012_k <- boot_net_sp_k.type.s(net_2012, 9, 40) %>%
  mutate(year = factor(rep("2012", n())))
spmet_2011_k <- boot_net_sp_k.type.s(net_2011, 4, 40) %>%
  mutate(year = factor(rep("2011", n())))  
spmet_2010_k <- boot_net_sp_k.type.s(net_2010, 3, 40) %>%
  mutate(year = factor(rep("2010", n())))

spmet_k <- bind_rows(spmet_2012_k, spmet_2011_k, spmet_2010_k)

spmet_mean_k <- spmet_k %>%
  group_by(site, year, elev.mean, metric) %>%
  summarize(value = mean(value))

### Network- and group-level
netmet_boot_2012_k <- boot_net_nt_k.type.s(net_2012, n = 9, iter = 40) %>%
  mutate(year = factor(rep("2012", n())))
netmet_boot_2011_k <- boot_net_nt_k.type.s(net_2011, n = 4, iter = 40)  %>%
  mutate(year = factor(rep("2011", n())))
netmet_boot_2010_k <- boot_net_nt_k.type.s(net_2010, n = 3, iter = 40) %>%
  mutate(year = factor(rep("2010", n())))

netmet_k <- bind_rows(netmet_boot_2012, 
                      netmet_boot_2011, 
                      netmet_boot_2010) %>%
  mutate(metric = factor(metric, levels = c("H2",
                                            "weighted.NODF",
                                            "niche.overlap.HL", 
                                            "C.score.HL",
                                            "generality.HL",
                                            "web.asymmetry",
                                            "weighted.connectance")))
```

# Plot k.type network analysis
```{r}
#### Per species values
ggplot(filter(spmet_2012_k,
              sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf")), 
       aes(elev.mean, value)) +
  geom_point() +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth() +
  facet_grid(metric~sp, scales = "free")

ggplot(filter(spmet_2011_k,
              sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf")), 
       aes(elev.mean, value)) +
  geom_point() +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth() +
  facet_grid(metric~sp, scales = "free")

ggplot(filter(spmet_2010_k,
              sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf")), 
       aes(elev.mean, value)) +
  geom_point() +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth() +
  facet_grid(metric~sp, scales = "free")

#### Community mean d
ggplot(filter(spmet_mean_k, metric == "d"), aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10), color = "black", linetype = "dashed", size = 0.75) +
  xlab("Elevation") +
  ylab("Community mean d'") +
  theme_light(16)

ggplot(filter(spmet_mean_k, metric == "d"), aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10), color = "black", linetype = "dashed", size = 0.75) +
  xlab("Elevation") +
  ylab("Community mean d'") +
  facet_wrap(~year) +
  theme_light(16)

gam_d <- gam(value ~ year + s(elev.mean, k = 10), 
             data = filter(spmet_k, metric == "d"),
             method = "REML")

gam_d.vis <- getViz(gam_d)
check.gamViz(gam_d.vis)
summary(gam_d.vis)

plot(gam_d.vis, allTerms = TRUE) +
  xlab("Elevation") 

#### Network level
ggplot(netmet_k, aes(elev.mean, value)) +
  geom_point() +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth() +
  facet_grid(metric~year, scales = "free")

##### H2
ggplot(filter(netmet_k, metric == "H2"), aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", linetype = "dashed", size = 0.75) +
  xlab("Elevation") +
  ylab("H2'") +
  theme_light(16)

gam_H2 <- gam(value ~ year + s(elev.mean), 
              data = filter(netmet_k, metric == "H2"),
              method = "REML")

gam_H2.vis <- getViz(gam_H2)
check.gamViz(gam_H2.vis)
summary(gam_H2.vis)

plot(gam_H2.vis, allTerms = TRUE) +
  xlab("Elevation") 

##### niche overlap
ggplot(filter(netmet_k, metric == "niche.overlap.HL"), aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", linetype = "dashed", size = 0.75) +
  #geom_smooth() +
  xlab("Elevation") +
  ylab("Niche overlap") +
  theme_light(16)

ggplot(filter(netmet_k, metric == "niche.overlap.HL"), aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", linetype = "dashed", size = 0.75) +
  #geom_smooth() +
  xlab("Elevation") +
  ylab("Niche overlap") +
  facet_wrap(~year) +
  theme_light(16)

gam_niche.overlap <- gam(value ~ year + s(elev.mean), 
                         data = filter(netmet_k, metric == "niche.overlap.HL"),
                         method = "REML")

gam_niche.overlap.vis <- getViz(gam_niche.overlap)
check.gamViz(gam_niche.overlap.vis)
summary(gam_niche.overlap.vis)

plot(gam_niche.overlap.vis, allTerms = TRUE) +
  xlab("Elevation") 


##### C.score
ggplot(filter(netmet_k, metric == "C.score.HL"), aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", linetype = "dashed", size = 0.75) +
  #geom_smooth() +
  xlab("Elevation") +
  ylab("C score") +
  theme_light(16)

ggplot(filter(netmet_k, metric == "C.score.HL"), aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", linetype = "dashed", size = 0.75) +
  #geom_smooth() +
  xlab("Elevation") +
  ylab("C score") +
  facet_wrap(~year) +
  theme_light(16)

gam_cscore <- gam(value ~ year + s(elev.mean), 
                  data = filter(netmet_k, metric == "C.score.HL"),
                  method = "REML")

gam_cscore.vis <- getViz(gam_cscore)
check.gamViz(gam_cscore.vis)
summary(gam_cscore.vis)

plot(gam_cscore.vis, allTerms = TRUE) +
  xlab("Elevation")
```


### Do BB species differ in their k.type associations?
```{r}
web_by_site <- net %>%
  group_by(site, bb.sp, k.type.s) %>% # group by site, date, bee species, plant species
  summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
  na.omit() %>%
  group_by(site) %>%
  filter(length(unique(bb.sp)) > 1 & length(unique(k.type.s)) > 1) %>%
  ungroup() %>%
  dplyr::select(higher = bb.sp, lower = k.type.s, webID = site, freq) %>% # rename columns to match bipartite's expectations
  data.frame() %>% # convert to data frame
  frame2webs()

web <- net %>%
  group_by(year, elev.class2, bb.sp, k.type.s) %>% # group by site, date, bee species, plant species
  summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
  na.omit() %>%
  unite(webID, c(year, elev.class2), sep = "_") %>%
  dplyr::select(higher = bb.sp, lower = k.type.s, webID, freq) %>% # rename columns to match bipartite's expectations
  data.frame() %>% # convert to data frame
  frame2webs()

# 2010 high
visweb(web[[1]])
plotweb(web[[1]])

# 2010 low
visweb(web[[2]])
plotweb(web[[2]])

# 2010 med
visweb(web[[3]])
plotweb(web[[3]])

# 2011 high
visweb(web[[4]])
plotweb(web[[4]])

# 2011 low
visweb(web[[5]])
plotweb(web[[5]])

# 2011 med
visweb(web[[6]])
plotweb(web[[6]])

# 2012 high
visweb(web[[7]])
plotweb(web[[7]])

# 2012 low
visweb(web[[8]])
plotweb(web[[8]])

# 2012 med
visweb(web[[9]])
plotweb(web[[9]])

```

### Does tongue length class predict k.type association? Remember that B. wurflenii is primarily a nectar robber.
```{r}
kugler_key <- read_delim("../processed_data/kugler_key.tsv", delim = "\t") %>%
  dplyr::select(1,2)

ggplot(bb_traits, aes(pbl.w)) +
  geom_histogram()

net_cont <- net %>%
  group_by(year, elev.class2, pbl.w.class, k.type.ss) %>%
  summarize(freq = n()) %>%
  group_by(year, elev.class2, pbl.w.class) %>%
  mutate(prop = freq/sum(freq))

net_cont2 <- net %>%
  group_by(year, elev.class2, bb.sp, k.type.ss) %>%
  summarize(freq = n()) %>%
  group_by(year, elev.class2, bb.sp) %>%
  mutate(prop = freq/sum(freq))

ggplot(net_cont, aes(pbl.w.class, freq, fill = k.type.ss)) +
  geom_col() +
  facet_grid(elev.class2 ~ year)

ggplot(net_cont, aes(pbl.w.class, prop, fill = k.type.ss)) +
  geom_col() +
  facet_grid(elev.class2 ~ year)

ggplot(net_cont2, aes(bb.sp, prop, fill = k.type.ss)) +
  geom_col() +
  facet_grid(elev.class2 ~ year)


chisq_mod <- chisq.test(net$pbl.w.class, net$k.type.ss)
chisq_mod

ggplot(net, aes(log(pbl.w), corolla.depth)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~year)

ggplot(net, aes(dayofyear, corolla.depth)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~year)

ggplot(net, aes(elev.mean, corolla.depth)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~year)

ggplot(net, aes(factor(pbl.w.class), corolla.depth)) +
  geom_boxplot() 
```

######## Old analyses that I probably don't care to use but don't want to delete yet #########

### Per-transect network analysis -- is it worth it? I would have to filter out the transects for which the networks are too small.
#### Network level
```{r eval=FALSE, include=FALSE}
per.transect_netmet <- net %>%
  group_by(site, dayofyear, gdd.cum, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
  summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
  group_by(site, dayofyear, gdd.cum, year) %>%
  mutate(size = sum(freq),
         bb.rich = length(unique(bb.sp)),
         fl.rich = length(unique(plant.sp)),
         total.rich = bb.rich + fl.rich) %>%
  filter(bb.rich > 3 & fl.rich > 3) %>%
  unite(webID, c(site, dayofyear, gdd.cum, year, size), sep = "_") %>%
  dplyr::select(higher = bb.sp, lower = plant.sp, webID, freq) %>% # rename columns to match bipartite's expectations
  data.frame() %>% # convert to data frame
  frame2webs() %>% # convert to to bipartite web
  map(networklevel, weighted = TRUE, level = "higher", index = c("niche overlap", "C score", "H2")) %>% # We probably care mainly about the BB level of the network
  data.frame() %>%
  t() %>%
  data.frame() %>%
  rownames_to_column() %>%
  dplyr::select(id = rowname, everything()) %>%
  as_tibble() %>%
  gather(metric, value, -id) %>%
  separate(id, into = c("site", "yday", "gdd.cum", "year", "size"), sep = "_") %>%
  mutate(size = as.numeric(size),
         year = factor(year),
         gdd.cum = as.numeric(gdd.cum),
         yday = as.numeric(yday),
         site = as.factor(site)) %>%
  left_join(site_data) %>%
  left_join(abund_div)

### Can we explain network metrics with abundance and diversity data?
# BB rich
ggplot(per.transect_netmet, aes(bb.rich, value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(metric ~ year, scales = "free")
# FL rich
ggplot(per.transect_netmet, aes(fl.rich, value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(metric ~ year, scales = "free")
# BB abund
ggplot(per.transect_netmet, aes(log(bb.abund), value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(metric ~ year, scales = "free")
# FL abund
ggplot(per.transect_netmet, aes(log(fl.abund), value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(metric ~ year, scales = "free")
# Per-capita FL abund
ggplot(per.transect_netmet, aes(log(fl.per.bb), value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(metric ~ year, scales = "free")


### Can we explain network metrics with elevation and/or time?
ggplot(per.transect_netmet, aes(elev.mean, value, color = log(size))) +
  geom_point() +
  geom_smooth() +
  facet_grid(metric ~ year, scales = "free")

ggplot(per.transect_netmet, aes(log(size), value, color = log(size))) +
  geom_point() +
  geom_smooth() +
  facet_grid(metric ~ year, scales = "free")

ggplot(filter(per.transect_netmet, size > 1), 
       aes(yday, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(metric ~ year, scales = "free")

ggplot(filter(per.transect_netmet, size > 10), 
       aes(gdd.cum, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(metric ~ year, scales = "free")

ggplot(filter(per.transect_netmet, size > 10), 
       aes(yday, elev.mean)) +
  geom_point(aes(color = value)) +
  geom_smooth() +
  facet_grid(metric ~ year, scales = "free")
```

#### Species level
```{r eval=FALSE, include=FALSE}
per.transect_spmet <- net %>%
  group_by(site, dayofyear, gdd.cum, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
  summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
  group_by(site, dayofyear, gdd.cum, year) %>%
  mutate(size = sum(freq),
         bb.rich = length(unique(bb.sp)),
         fl.rich = length(unique(plant.sp)),
         total.rich = bb.rich + fl.rich) %>%
  filter(bb.rich > 5 & fl.rich > 5) %>%
  unite(webID, c(site, dayofyear, gdd.cum, year, size), sep = "_") %>%
  dplyr::select(higher = bb.sp, lower = plant.sp, webID, freq) %>% # rename columns to match bipartite's expectations
  data.frame() %>% # convert to data frame
  frame2webs() %>% # convert to to bipartite web
  map(specieslevel, level = "higher", index = "d") %>% # We probably care mainly about the BB level of the network
  map(rownames_to_column) %>%
  map(as_tibble) %>%
  bind_rows(.id = "id") %>% # Collapse list into big data frame
  rename(sp = rowname) %>%
  gather(metric, value, -c(sp, id)) %>% # get mean value across iterations for each metric
  separate(id, sep = "_", into = c("site", "yday", "gdd.cum", "year", "size")) %>%
  mutate(size = as.numeric(size),
         year = factor(year),
         gdd.cum = as.numeric(gdd.cum),
         yday = as.numeric(yday)) %>%
  left_join(site_data) %>%
  left_join(bb_traits, by = c("sp" = "bb.sp")) %>%
  left_join(abund_div)

per.transect_spmet.mean <- per.transect_spmet %>%
  group_by(site, yday, gdd.cum, year) %>%
  summarize(value = mean(value)) %>%
  left_join(site_data)


ggplot(filter(per.transect_spmet), 
       aes(log(fl.abund), value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(. ~ year, scales = "free")

ggplot(filter(per.transect_spmet), 
       aes(elev.mean, size)) +
  geom_point() +
  geom_smooth() +
  facet_grid(. ~ year, scales = "free")

ggplot(filter(per.transect_spmet.mean), 
       aes(yday, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(. ~ year, scales = "free")

ggplot(filter(per.transect_spmet), 
       aes(elev.mean, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(. ~ year, scales = "free")

ggplot(filter(per.transect_spmet, size >= 10), 
       aes(gdd.cum, value, color = log(size))) +
  geom_point() +
  geom_smooth() +
  facet_grid(metric ~ year, scales = "free")

ggplot(filter(per.transect_spmet, size >= 10), 
       aes(yday, value, color = log(size))) +
  geom_point() +
  geom_smooth() +
  facet_grid(metric ~ year, scales = "free")
```


### Phase x elev network analysis -- this only makes sense in 2012, where we have sufficient sampling across the whole colony lifecycle
```{r eval=FALSE, include=FALSE}
### Prepare pxe abundance and diversity table
pxe_abund_div <- abund_div %>%
  dplyr::select(-gdd.cum, -yday, -week, -month, -date, -elev.class2) %>%
  gather(ad.metric, ad.value, -c(site, year, phase, elev.mean)) %>%
  group_by(site, phase, year, ad.metric) %>%
  summarize(ad.value = mean(ad.value, na.rm = TRUE))


pxe_net <- boot_net_nt2(filter(net, year == 2012), 3, 20) %>%
  full_join((pxe_abund_div))


ggplot(filter(pxe_net, metric == "H2"), aes(elev.mean, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(year ~ phase, scales = "free") +
  xlab("Elevation") +
  ylab("H2")

ggplot(filter(pxe_net, metric == "C.score.HL"), aes(elev.mean, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(year ~ phase, scales = "free") +
  xlab("Elevation") +
  ylab("C.score")

ggplot(filter(pxe_net, metric == "niche.overlap.HL"), aes(elev.mean, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(year ~ phase, scales = "free") +
  xlab("Elevation") +
  ylab("Niche overlap")


pxe_sp <- boot_net_sp2(net, 4, 20)  %>%
  full_join((pxe_abund_div))

ggplot(filter(pxe_sp, 
              metric == "d" &
                sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf", "psit") &
                ad.metric == "fl.rich"), 
       aes(ad.value, value, color = year)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(sp ~ phase, scales = "free") +
  ylim(c(0,1)) +
  xlab("Floral richness") +
  ylab("d'")

pxe_sp_mean <- pxe_sp %>%
  filter(sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf", "psit")) %>%
  group_by(site, year, phase, elev.mean, metric) %>%
  summarize(value = mean(value, na.rm = TRUE))

ggplot(filter(pxe_sp_mean, metric == "d"), aes(elev.mean, value)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  facet_grid(year ~ phase, scales = "free") +
  xlab("Elevation") +
  ylab("Mean d'")

ggplot(filter(pxe_sp_mean, metric == "species.strength"), aes(elev.mean, value)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  facet_grid(year ~ phase, scales = "free") +
  xlab("Elevation") +
  ylab("Mean species strength'")

ggplot(filter(pxe_sp, 
              metric == "species.strength",
              sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf", "psit")), 
       aes(elev.mean, value, color = year)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  facet_grid(sp ~ phase, scales = "free") +
  xlab("Elevation") +
  ylab("Species strength")

ggplot(filter(pxe_sp, 
              metric == "d",
              sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf", "psit")), 
       aes(elev.mean, value, color = year)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  facet_grid(sp ~ phase, scales = "free") +
  ylim(c(0,1)) +
  xlab("Elevation") +
  ylab("d'")

ggplot(filter(pxe_sp, 
              metric == "d",
              sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf", "psit")), 
       aes(elev.mean, value)) +
  geom_point(aes(color = year)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  facet_grid(sp ~ phase, scales = "free") +
  ylim(c(0,1)) +
  xlab("Elevation") +
  ylab("d'")

ggplot(filter(pxe_sp, 
              metric == "d",
              sp %in% c("bss", "hort", "pasc", "prat", "soro", "wurf", "psit")), 
       aes(elev.mean, value, color = phase)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
  facet_grid(year ~ sp, scales = "free") +
  theme_light() +
  xlab("Elevation") +
  ylab("d'")
```

```{r}

# network level PxE
boot_net_nt_pxe <- function(net, n, iter) {
  
  out <- data.frame()
  
  t_filter <- net %>%
    dplyr::select(site, date, phase, year) %>%
    unique() %>%
    group_by(site, phase, year) %>%
    mutate(transects = n()) %>%
    filter(transects >= n)
  
  for(i in 1:iter) {
    temp <- net %>%
      semi_join(t_filter) %>%
      dplyr::select(site, date, phase, year) %>%
      unique() %>%
      group_by(site, phase, year) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, phase, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
      summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
      group_by(site, phase, year) %>%
      filter(length(unique(bb.sp)) > 1 & length(unique(plant.sp)) > 1) %>%
      ungroup() %>%
      unite(webID, c(site, phase, year), sep = "_") %>%
      dplyr::select(higher = bb.sp, lower = plant.sp, webID, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      map(networklevel, weighted = TRUE, level = "higher",
          index = c("web asymmetry", "weighted NODF", "H2", "C score", "niche overlap", "generality")) %>% 
      data.frame() %>%
      t() %>%
      data.frame() %>%
      rownames_to_column() %>%
      dplyr::select(id = rowname, everything()) %>%
      as_tibble() %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out  %>%
    gather(metric, value, -c(id, iteration)) %>% # get mean value across iterations for each metric
    group_by(id, metric) %>%
    summarize(value = mean(value)) %>%
    separate(id, sep = "_", into = c("site", "phase", "year")) %>%
    mutate(phase = factor(phase, levels = c("founding", "buildup", "reproductive"))) %>%
    left_join(site_data)
}

# species level version 2
boot_net_sp2 <- function(net, n, iter) {
  
  out <- data.frame()
  
  t_filter <- net %>%
    dplyr::select(site, date, phase, year) %>%
    unique() %>%
    group_by(site, phase, year) %>%
    mutate(transects = n()) %>%
    filter(transects >= n)
  
  for(i in 1:iter) {
    temp <- net %>%
      semi_join(t_filter) %>%
      dplyr::select(site, date, phase, year) %>%
      unique() %>%
      group_by(site, phase, year) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, phase, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
      summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
      group_by(site, phase, year) %>%
      filter(length(unique(bb.sp)) > 1 & length(unique(plant.sp)) > 1) %>%
      ungroup() %>%
      unite(webID, c(site, phase, year), sep = "_") %>%
      dplyr::select(higher = bb.sp, lower = plant.sp, webID, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      map(specieslevel, level = "higher", index = c("d", "species strength")) %>%
      map(rownames_to_column) %>%
      map(as_tibble) %>%
      bind_rows(.id = "id") %>% # Collapse list into big data frame
      rename(sp = rowname) %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out %>%
    gather(metric, value, -c(sp, id, iteration)) %>% # get mean value across iterations for each metric
    group_by(sp, id, metric) %>%
    summarize(value = mean(value)) %>%
    separate(id, sep = "_", into = c("site", "phase", "year")) %>%
    mutate(phase = factor(phase, levels = c("founding", "buildup", "reproductive"))) %>%
    left_join(site_data) %>%
    left_join(bb_traits, by = c("sp" = "bb.sp"))
}

```


#### plants
spmet_2012_FL <- boot_net_sp_FL(net_2012, 9, 40) %>%
  mutate(year = factor(rep("2012", n())))
spmet_2011_FL <- boot_net_sp_FL(net_2011, 4, 40) %>%
  mutate(year = factor(rep("2011", n())))  
spmet_2010_FL <- boot_net_sp_FL(net_2010, 3, 40) %>%
  mutate(year = factor(rep("2010", n())))

spmet_FL <- bind_rows(spmet_2012_FL, spmet_2011_FL, spmet_2010_FL)

spmet_mean_FL <- spmet_FL %>%
  group_by(site, year, elev.mean, metric) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(year = factor(year))


### Demonstrate that the distribution of interaction frequencies does not differ across elevation classes; in other words, no systematic bias is introduced to the network analysis due to uneven frequency of species with few observations
net_freq <- net %>%
  group_by(site, dayofyear, year, bb.sp, elev.mean, elev.class2) %>%
  summarize(freq = n())

ggplot(net_freq, aes(elev.class2, log(freq))) +
  geom_boxplot() +
  facet_wrap(~year)


### Let's try to make a bootstrapped version of this to deal with the unevenness in # transects
```{r}
betalinkr_boot <- function(net, n, iter) {
  out <- data.frame()
  
  for(i in 1:iter) {
    
    temp <- net %>%
      dplyr::select(site, date, year) %>%
      unique() %>%
      group_by(site, year) %>%
      sample_n(n) %>%
      left_join(net) %>%
      group_by(site, bb.sp, plant.sp) %>%
      summarize(freq = n()) %>% 
      dplyr::select(higher = bb.sp, lower = plant.sp, webID = site, freq) %>% # rename columns to match bipartite's expectations
      data.frame() %>% # convert to data frame
      frame2webs() %>% # convert to to bipartite web
      webs2array() %>% # convert to array
      betalinkr_multi(partition.st = TRUE, partition.rr = TRUE) %>%
      as_tibble() %>%
      rename(site1 = i) %>%
      rename(site2 = j) %>%
      mutate(iteration = rep(i, n()))
    
    out <- bind_rows(out, temp)
  } 
  out <- out %>%
    gather(metric, value, -c(site1, site2, iteration)) %>% # get mean value across iterations for each metric
    group_by(site1, site2, metric) %>%
    summarize(value = mean(value))
}
```

### Let's take betalinkr_boot for a spin
```{r}
beta1 <- betalinkr_boot(net, 3, 10)
```

netmet_sitedate <- site_date_net(net, min_higher = 2, min_lower = 2) %>%
  mutate(site = factor(site),
         date = ymd(date),
         yday = yday(date),
         year = factor(year(date)))

ggplot(netmet_sitedate, aes(elev.mean, niche.overlap.HL)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~year)

ggplot(netmet_sitedate, aes(elev.mean, niche.overlap.HL)) +
  geom_point() +
  geom_smooth()

ggplot(netmet_sitedate, aes(yday, niche.overlap.HL, color = year)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
  facet_wrap(~site, scales = "free")

bigbam_2010 <- bamV(niche.overlap.HL ~ s(site, bs = "re") +
                      te(yday, elev.mean, 
                         bs= c("tp", "tp"), 
                         k=c(8, 8), m=1),
                    data = filter(netmet_sitedate, year == 2010))

check.gamViz(bigbam_2010)
bigbam_2010.sum <- summary(bigbam_2010)
print(plot(bigbam_2010, allTerms = TRUE), pages = 1)

bigbam_2011 <- bamV(niche.overlap.HL ~ s(site, bs = "re") +
                      te(yday, elev.mean, 
                         bs= c("tp", "tp"), 
                         k=c(8, 8), m=1),
                    data = filter(netmet_sitedate, year == 2011))

check.gamViz(bigbam_2011)
bigbam_2011.sum <- summary(bigbam_2011)
print(plot(bigbam_2011, allTerms = TRUE), pages = 1)

bigbam_2012 <- bamV(niche.overlap.HL ~ s(site, bs = "re") +
                      te(yday, elev.mean, 
                         bs= c("tp", "tp"), 
                         k=c(8, 8), m=1),
                    data = filter(netmet_sitedate, year == 2012))

check.gamViz(bigbam_2012)
bigbam_2012.sum <- summary(bigbam_2012)
print(plot(bigbam_2012, allTerms = TRUE), pages = 1)

### Transect-level analysis
```{r}
beta_site.date.2010 <- betalinkr_multi(array_site.date.2010, 
                                       partition.st = TRUE, 
                                       partition.rr = TRUE) %>%
  as_tibble() %>%
  separate(i, into = c("site1", "year1", "yday1"), sep = "_") %>%
  separate(j, into = c("site2", "year2", "yday2"), sep = "_") %>%
  left_join(site_data, by = c("site1" = "site"))%>%
  dplyr::select(site1, year1, yday1, site2, year2, yday2, S, OS, WN, ST, 
                ST.l, ST.h, ST.lh, WN.repl, OS.repl, WN.rich, OS.rich, 
                elev1 = elev.mean, man1 = management, lon1 = lon, lat1 = lat) %>%
  left_join(site_data, by = c("site2" = "site")) %>%
  dplyr::select(site1, year1, yday1, site2, year2, yday2, S, OS, WN, ST, 
                ST.l, ST.h, ST.lh, WN.repl, OS.repl, WN.rich, OS.rich, 
                elev1, elev2 = elev.mean, man1, man2 = management, lon1, lat1, lon2 = lon, lat2 = lat) %>%
  mutate(site1 = factor(site1),
         site2 = factor(site2),
         year1 = factor(year1),
         year2 = factor(year2),
         yday1 = as.numeric(yday1),
         yday2 = as.numeric(yday2),
         elev.diff = abs(elev1 - elev2),
         mean.elev = (elev1 + elev2)/2,
         yday.diff = abs(yday1 - yday2),
         man.diff = case_when(
           man1 == man2 ~ FALSE,
           man1 != man2 ~ TRUE
         )) %>%
  gather(metric, value, -c(site1, year1, yday1, site2, year2, yday2,
                           elev1, elev2, elev.diff, yday.diff, 
                           man1, man2, man.diff, mean.elev, 
                           lon1, lat1, lon2, lat2)) %>%
  mutate(metric.class = case_when(
    metric %in% c("WN", "ST", "OS", "S") ~ "aggregate",
    metric %in% c("WN.repl", "WN.rich") ~ "WN.partition",
    metric %in% c("ST.h", "ST.l", "ST.lh") ~ "ST.partition",
    metric %in% c("OS.repl", "OS.rich") ~ "OS.partition"
  ))

beta_site.date.2011 <- betalinkr_multi(array_site.date.2011, 
                                       partition.st = TRUE, 
                                       partition.rr = TRUE) %>%
  as_tibble() %>%
  separate(i, into = c("site1", "year1", "yday1"), sep = "_") %>%
  separate(j, into = c("site2", "year2", "yday2"), sep = "_") %>%
  left_join(site_data, by = c("site1" = "site"))%>%
  dplyr::select(site1, year1, yday1, site2, year2, yday2, S, OS, WN, ST, 
                ST.l, ST.h, ST.lh, WN.repl, OS.repl, WN.rich, OS.rich, 
                elev1 = elev.mean, man1 = management, lon1 = lon, lat1 = lat) %>%
  left_join(site_data, by = c("site2" = "site")) %>%
  dplyr::select(site1, year1, yday1, site2, year2, yday2, S, OS, WN, ST, 
                ST.l, ST.h, ST.lh, WN.repl, OS.repl, WN.rich, OS.rich, 
                elev1, elev2 = elev.mean, man1, man2 = management, lon1, lat1, lon2 = lon, lat2 = lat) %>%
  mutate(site1 = factor(site1),
         site2 = factor(site2),
         year1 = factor(year1),
         year2 = factor(year2),
         yday1 = as.numeric(yday1),
         yday2 = as.numeric(yday2),
         elev.diff = abs(elev1 - elev2),
         mean.elev = (elev1 + elev2)/2,
         yday.diff = abs(yday1 - yday2),
         man.diff = case_when(
           man1 == man2 ~ FALSE,
           man1 != man2 ~ TRUE
         )) %>%
  gather(metric, value, -c(site1, year1, yday1, site2, year2, yday2,
                           elev1, elev2, elev.diff, yday.diff, 
                           man1, man2, man.diff, mean.elev, 
                           lon1, lat1, lon2, lat2)) %>%
  mutate(metric.class = case_when(
    metric %in% c("WN", "ST", "OS", "S") ~ "aggregate",
    metric %in% c("WN.repl", "WN.rich") ~ "WN.partition",
    metric %in% c("ST.h", "ST.l", "ST.lh") ~ "ST.partition",
    metric %in% c("OS.repl", "OS.rich") ~ "OS.partition"
  ))

beta_site.date.2012 <- betalinkr_multi(array_site.date.2012, 
                                       partition.st = TRUE, 
                                       partition.rr = TRUE) %>%
  as_tibble() %>%
  separate(i, into = c("site1", "year1", "yday1"), sep = "_") %>%
  separate(j, into = c("site2", "year2", "yday2"), sep = "_") %>%
  left_join(site_data, by = c("site1" = "site"))%>%
  dplyr::select(site1, year1, yday1, site2, year2, yday2, S, OS, WN, ST, 
                ST.l, ST.h, ST.lh, WN.repl, OS.repl, WN.rich, OS.rich, 
                elev1 = elev.mean, man1 = management, lon1 = lon, lat1 = lat) %>%
  left_join(site_data, by = c("site2" = "site")) %>%
  dplyr::select(site1, year1, yday1, site2, year2, yday2, S, OS, WN, ST, 
                ST.l, ST.h, ST.lh, WN.repl, OS.repl, WN.rich, OS.rich, 
                elev1, elev2 = elev.mean, man1, man2 = management, lon1, lat1, lon2 = lon, lat2 = lat) %>%
  mutate(site1 = factor(site1),
         site2 = factor(site2),
         year1 = factor(year1),
         year2 = factor(year2),
         yday1 = as.numeric(yday1),
         yday2 = as.numeric(yday2),
         elev.diff = abs(elev1 - elev2),
         mean.elev = (elev1 + elev2)/2,
         yday.diff = abs(yday1 - yday2),
         man.diff = case_when(
           man1 == man2 ~ FALSE,
           man1 != man2 ~ TRUE
         )) %>%
  gather(metric, value, -c(site1, year1, yday1, site2, year2, yday2,
                           elev1, elev2, elev.diff, yday.diff, 
                           man1, man2, man.diff, mean.elev, 
                           lon1, lat1, lon2, lat2)) %>%
  mutate(metric.class = case_when(
    metric %in% c("WN", "ST", "OS", "S") ~ "aggregate",
    metric %in% c("WN.repl", "WN.rich") ~ "WN.partition",
    metric %in% c("ST.h", "ST.l", "ST.lh") ~ "ST.partition",
    metric %in% c("OS.repl", "OS.rich") ~ "OS.partition"
  ))

beta_site.date <- bind_rows(beta_site.date.2010, beta_site.date.2011, beta_site.date.2012)
beta_site.date_gdm.prep <- beta_site.date %>%
  filter(metric == "WN")

write_csv(beta_site.date, "../analyses/output/beta_site.date.csv")
beta_site.date <- read_csv("../analyses/output/beta_site.date.csv")
```


# Beta diversity by elevation difference
```{r}

### Okay, we'll just do partial mantel tests to partial out the effect of transects.diff; and we'll only do that for the aggregate measures because it seems a bit silly to doi it for the partitions  
WN_diff <- beta_site %>%
  filter(metric == "WN") %>%
  dplyr::select(site1, site2, value) %>%
  spread(site2, value) %>%
  dist()

WN_diff <- beta_site %>%
  filter(metric == "WN") %>%
  dplyr::select(site1, site2, value) %>%
  spread(site2, value) %>%
  dist()

ST_diff <- beta_site %>%
  filter(metric == "ST") %>%
  dplyr::select(site1, site2, value) %>%
  spread(site2, value) %>%
  dist()

OS_diff <- beta_site %>%
  filter(metric == "OS") %>%
  dplyr::select(site1, site2, value) %>%
  spread(site2, value) %>%
  dist()

elev_diff <- beta_site %>%
  dplyr::select(site1, site2, elev.diff) %>%
  unique() %>%
  spread(site2, elev.diff) %>%
  dist()

transects_diff <- beta_site %>%
  dplyr::select(site1, site2, transects.diff) %>%
  unique() %>%
  spread(site2, transects.diff) %>%
  dist()

WN_mantel <- mantel.partial(elev_diff, WN_diff, transects_diff)
ST_mantel <- mantel.partial(elev_diff, ST_diff, transects_diff)
OS_mantel <- mantel.partial(elev_diff, OS_diff, transects_diff)

### 2011
#### Let's explore whether the different transformation options I'm weighing have any effect on the overall pattern that we end up visualizing

#### Option A: log(n + 0.001) with gaussian
weighted_fl_abundance.2011.1 <- bam(log(fl.abund.wt + 0.01) ~
                                      s(bb.sp, bs = "re") +
                                      #s(management, bs = "re") +
                                      #s(site, bs = "re") +
                                      te(yday, elev.mean, 
                                         by = bb.sp,
                                         #id = 0,
                                         bs= c("gp", "tp"), 
                                         k=c(8, 8), 
                                         m=2),
                                    data = filter(fl_abund_wt, year == 2011),
                                    discrete = TRUE,
                                    family = "gaussian",
                                    select = TRUE,
                                    method = "fREML") %>% getViz(nsim = 500)

check.gamViz(weighted_fl_abundance.2011.1)
#weighted_fl_abundance.2011.sum.1 <- summary(weighted_fl_abundance.2011.1)
#print(plot(weighted_fl_abundance.2011.1, allTerms = TRUE), pages = 4)
weighted_fl_abundance.2011.ck.1 <- check2D(weighted_fl_abundance.2011.1, 
                                           x1 = "yday", x2 = "elev.mean")

#### Option B: sqrt(n) with tweedie
weighted_fl_abundance.2011.2 <- bam(sqrt(fl.abund.wt) ~
                                      s(bb.sp, bs = "re") +
                                      #s(management, bs = "re") +
                                      #s(site, bs = "re") +
                                      te(yday, elev.mean, 
                                         by = bb.sp,
                                         #id = 0,
                                         bs= c("gp", "tp"), 
                                         k=c(8, 8), 
                                         m=2),
                                    data = filter(fl_abund_wt, year == 2011),
                                    discrete = TRUE,
                                    family = "tw",
                                    select = TRUE,
                                    method = "fREML") %>% getViz(nsim = 500)

check.gamViz(weighted_fl_abundance.2011.2)
#weighted_fl_abundance.2011.sum.2 <- summary(weighted_fl_abundance.2011.2)
#print(plot(weighted_fl_abundance.2011.2, allTerms = TRUE), pages = 4)
weighted_fl_abundance.2011.ck.2 <- check2D(weighted_fl_abundance.2011.2, 
                                           x1 = "yday", x2 = "elev.mean")


#### Option C: sqrt(sqrt.n) with tweedie
weighted_fl_abundance.2011.3 <- bam(sqrt(fl.abund.wt.sqrt) ~
                                      s(bb.sp, bs = "re") +
                                      #s(management, bs = "re") +
                                      #s(site, bs = "re") +
                                      te(yday, elev.mean, 
                                         by = bb.sp,
                                         #id = 0,
                                         bs= c("gp", "tp"), 
                                         k=c(8, 8), 
                                         m=2),
                                    data = filter(fl_abund_wt, year == 2011),
                                    discrete = TRUE,
                                    family = "tw",
                                    select = TRUE,
                                    method = "fREML") %>% getViz(nsim = 500)

check.gamViz(weighted_fl_abundance.2011.3)
#weighted_fl_abundance.2011.sum.3 <- summary(weighted_fl_abundance.2011.3)
#print(plot(weighted_fl_abundance.2011, allTerms = TRUE), pages = 4)
weighted_fl_abundance.2011.ck.3 <- check2D(weighted_fl_abundance.2011.3, 
                                           x1 = "yday", x2 = "elev.mean")

#### Option D: sqrt(bin.n) with tweedie
weighted_fl_abundance.2011.4 <- bam(sqrt(fl.abund.wt.bin) ~
                                      s(bb.sp, bs = "re") +
                                      #s(management, bs = "re") +
                                      #s(site, bs = "re") +
                                      te(yday, elev.mean, 
                                         by = bb.sp,
                                         #id = 0,
                                         bs= c("gp", "tp"), 
                                         k=c(8, 8), 
                                         m=2),
                                    data = filter(fl_abund_wt, year == 2011),
                                    discrete = TRUE,
                                    family = "tw",
                                    select = TRUE,
                                    method = "fREML") %>% getViz(nsim = 500)

check.gamViz(weighted_fl_abundance.2011.4)
#weighted_fl_abundance.2011.sum.4 <- summary(weighted_fl_abundance.2011.4)
#print(plot(weighted_fl_abundance.2011, allTerms = TRUE), pages = 4)
weighted_fl_abundance.2011.ck.4 <- check2D(weighted_fl_abundance.2011.4, 
                                           x1 = "yday", x2 = "elev.mean")

#### Option E: log(sqrt.n + 0.001) with gaussian
weighted_fl_abundance.2011.5 <- bam(log(fl.abund.wt.sqrt + 0.001) ~
                                      s(bb.sp, bs = "re") +
                                      #s(management, bs = "re") +
                                      #s(site, bs = "re") +
                                      te(yday, elev.mean, 
                                         by = bb.sp,
                                         #id = 0,
                                         bs= c("gp", "tp"), 
                                         k=c(8, 8), 
                                         m=2),
                                    data = filter(fl_abund_wt, year == 2011),
                                    discrete = TRUE,
                                    family = "gaussian",
                                    select = TRUE,
                                    method = "fREML") %>% getViz(nsim = 500)

check.gamViz(weighted_fl_abundance.2011.5)
#weighted_fl_abundance.2011.sum.5 <- summary(weighted_fl_abundance.2011.5)
#print(plot(weighted_fl_abundance.2011, allTerms = TRUE), pages = 4)
weighted_fl_abundance.2011.ck.5 <- check2D(weighted_fl_abundance.2011.5, 
                                           x1 = "yday", x2 = "elev.mean")

#### Option F: log(bin.n + 0.01) with gaussian
weighted_fl_abundance.2011.6 <- bam(log(fl.abund.wt.bin + 0.01) ~
                                      s(bb.sp, bs = "re") +
                                      #s(management, bs = "re") +
                                      #s(site, bs = "re") +
                                      te(yday, elev.mean, 
                                         by = bb.sp,
                                         #id = 0,
                                         bs= c("gp", "tp"), 
                                         k=c(8, 8), 
                                         m=2),
                                    data = filter(fl_abund_wt, year == 2011),
                                    discrete = TRUE,
                                    family = "gaussian",
                                    select = TRUE,
                                    method = "fREML") %>% getViz(nsim = 500)

check.gamViz(weighted_fl_abundance.2011.6)
#weighted_fl_abundance.2011.sum.6 <- summary(weighted_fl_abundance.2011.6)
#print(plot(weighted_fl_abundance.2011.6, allTerms = TRUE), pages = 4)
weighted_fl_abundance.2011.ck.6 <- check2D(weighted_fl_abundance.2011.6, 
                                           x1 = "yday", x2 = "elev.mean")
```

### Since we clearly have structured collinearity in our preditors, i.e. elevation shapes abundance and diversity which shape network structure, perhaps SEM could be appropriate? The main question is whether it can handle the potentially complex (nonlinear) relationships between elevation and diversity/abundance, and whether it can handle the hierarchical structure of my data.
```{r}
# Create scales variables; this is important since SEM is based on covariance (right?)
# For starters, let's just try to model niche.overlap
netmet.scaled <- netmet %>%
  select(site, year, niche.overlap.HL, elev.mean, fl.rich, bb.rich, log.fl.abund, log.bb.abund) %>%
  mutate(elev.mean.2 = elev.mean^2,
         elev.mean.3 = elev.mean^3) %>%
  mutate(niche.overlap.HL = scale(niche.overlap.HL),
         elev.mean = scale(elev.mean),
         elev.mean.2 = scale(elev.mean),
         elev.mean.3 = scale(elev.mean^3),
         fl.rich = scale(fl.rich),
         log.fl.abund = scale(log.fl.abund),
         log.bb.abund = scale(log.bb.abund)) %>%
  as.matrix() %>%
  as_tibble() %>%
  mutate(niche.overlap.HL = as.numeric(niche.overlap.HL),
         elev.mean = as.numeric(elev.mean),
         elev.mean.2 = as.numeric(elev.mean^2),
         elev.mean.3 = as.numeric(elev.mean^3),
         fl.rich = as.numeric(fl.rich),
         bb.rich = as.numeric(bb.rich),
         log.fl.abund = as.numeric(log.fl.abund),
         log.bb.abund = as.numeric(log.bb.abund),
         site = factor(site),
         year = factor(year)) %>%
  na.omit()

# Pairs plot to determine best model form
predictors.scaled <- netmet.scaled %>%
  select(elev.mean, elev.mean.2, elev.mean.3, bb.rich, log.bb.abund, fl.rich, log.fl.abund, niche.overlap.HL)

ggpairs(predictors.scaled)

# Diversity and abundance metrics clearly have a quadratic relationship with elevation

# Call models
# To avoid figuring out how to handle year and site as random effects, let's start with separate year models
netmet.scaled.2011 <- filter(netmet.scaled, year == 2011)

mod_niche.overlap.2011 <- psem(
  lm(niche.overlap.HL ~ elev.mean + elev.mean.2 + elev.mean.3 + fl.rich + log.bb.abund, 
     netmet.scaled.2011),
  lm(fl.rich ~ elev.mean + elev.mean.2, 
     netmet.scaled.2011),
  lm(log.bb.abund ~ elev.mean + elev.mean.2 + fl.rich, 
     netmet.scaled.2011)
)

summary(mod_niche.overlap.2011)
plot(mod_niche.overlap.2011)
```

### But let's zoom in a bit closer on d' and see if we can model it per-species
```{r}
gam_d.sp_GI <- gam(d ~ year + 
                     #s(site, bs = "re") +
                     #management +
                     s(bb.sp, bs = "re") +
                     s(elev.mean, bs = "tp", k = 8) +
                     s(elev.mean, by = bb.sp, bs = "tp", k = 8) +
                     s(fl.rich, k = 8) +
                     s(log.fl.abund, k = 8) +
                     s(log.bb.sp.abund, k = 5),
                   select = TRUE,
                   data = spmet,
                   method = "REML") %>% getViz()

check.gamViz(gam_d.sp_GI)
concurvity(gam_d.sp_GI, full = TRUE)
concurvity(gam_d.sp_GI, full = FALSE)
summary(gam_d.sp_GI)
print(plot(gam_d.sp_GI, allTerms = TRUE), pages = 4)

gam_d.sp_G <- gam(d ~ year + 
                    s(site, bs = "re") +
                    management +
                    s(bb.sp, bs = "re") +
                    s(elev.mean, bs = "tp", k = 8) +
                    s(fl.rich, k = 8) +
                    s(log.fl.abund, k = 8) +
                    s(log.bb.sp.abund, k = 5),
                  select = TRUE,
                  data = spmet,
                  method = "REML") %>% getViz()

check.gamViz(gam_d.sp_G)
concurvity(gam_d.sp_G, full = TRUE)
concurvity(gam_d.sp_G, full = FALSE)
summary(gam_d.sp_G)
print(plot(gam_d.sp_G, allTerms = TRUE), pages = 1)
check_model(gam_d.sp_G)
modelG_elev.mean <- ggpredict(gam_d.sp_G, terms = "elev.mean")
plot(modelG_elev.mean)

# This simpler varying-smooth model explains just as much variance as the more complex GI model, so it should be preferred. 
gam_d.sp_GS <- gam(d ~ year + 
                     management +
                     s(site, bs = "re") +
                     s(elev.mean, bs = "tp", k = 8) +
                     s(elev.mean, bb.sp, bs = "fs", k = 8) +
                     s(log.fl.abund, k = 8) +
                     s(fl.rich, k = 8) +
                     s(log.bb.sp.abund, k = 5),
                   select = TRUE,
                   data = spmet,
                   method = "REML") %>% getViz()

check.gamViz(gam_d.sp_GS)
concurvity(gam_d.sp_GS, full = TRUE)
concurvity(gam_d.sp_GS, full = FALSE)
summary(gam_d.sp_GS)
print(plot(gam_d.sp_GS, allTerms = TRUE), pages = 1)
check_model(gam_d.sp_GS)
plot(ggpredict(gam_d.sp_GS))

# What about I? Does it do any better than S?
gam_d.sp_I <- gam(d ~ year + 
                    #s(site, bs = "re") +
                    #management +
                    #s(bb.transects) +
                    s(bb.sp, bs = "re") +
                    s(elev.mean, by = bb.sp, bs = "tp", k = 10),
                  #s(fl.rich, k = 8) +
                  #s(log.fl.abund, k = 8) +
                  #s(bb.sp.abund, k = 5),
                  select = TRUE,
                  data = filter(spmet, bb.sp %in% c("pasc", "prat", "bss", "soro", "wurf")),
                  method = "REML") %>% getViz()

check.gamViz(gam_d.sp_I)
concurvity(gam_d.sp_I, full = TRUE)
concurvity(gam_d.sp_I, full = FALSE)
summary(gam_d.sp_I)
print(plot(gam_d.sp_I, allTerms = TRUE), pages = 4)

AIC(gam_d.sp_I, gam_d.sp_S, gam_d.sp_G, gam_d.sp_GI, gam_d.sp_GS)

ggplot(spmet, aes(log.bb.sp.abund)) +
  geom_histogram() +
  theme_modern()

ggplot(spmet, aes(log.bb.sp.abund)) +
  geom_histogram() +
  theme_blackboard()

ggplot(spmet, aes(elev.mean, d)) +
  geom_point(aes(color = bb.sp)) +
  geom_smooth() +
  theme_blackboard() +
  facet_wrap(~year)

# We lose virtually no explanatory power in going from GI, to GS, to S. So, S is the preferred model right now.
gam_d.sp_S <- gam(d ~ year + 
                    s(site, bs = "re") +
                    management +
                    s(elev.mean, bb.sp, bs = "fs", k = 8) +
                    s(fl.rich, k = 8) +
                    s(log.fl.abund, k = 8) +
                    s(log.bb.sp.abund, k = 5),
                  select = TRUE,
                  data = spmet,
                  method = "REML") %>% getViz()

check.gamViz(gam_d.sp_S)
concurvity(gam_d.sp_S, full = TRUE)
concurvity(gam_d.sp_S, full = FALSE)
summary(gam_d.sp_S)
print(plot(gam_d.sp_S, allTerms = TRUE), pages = 1)


###  There remains a lot of noise in my multi-year models above, and a possible explanation for that is that the influence of year extends beyond just variation in the intercept. What if year affects the shape of Y ~ X? In these year subset models, I can drop both 'year' and 's(site, bs = "re")'.

gam_d.sp_GI.2010 <- gam(d ~ management +
                          s(bb.sp, bs = "re") +
                          s(elev.mean, bs = "tp", k = 8) +
                          s(elev.mean, by = bb.sp, bs = "tp", k = 8) +
                          s(fl.rich, k = 8) +
                          s(log.fl.abund, k = 8) +
                          s(log.bb.sp.abund, k = 5),
                        select = TRUE,
                        data = filter(spmet, year == 2010),
                        method = "REML") %>% getViz()

check.gamViz(gam_d.sp_GI.2010)
concurvity(gam_d.sp_GI.2010, full = TRUE)
concurvity(gam_d.sp_GI.2010, full = FALSE)
summary(gam_d.sp_GI.2010)
print(plot(gam_d.sp_GI.2010, allTerms = TRUE), pages = 4)

gam_d.sp_GI.2011 <- gam(d ~ 
                          #management +
                          s(bb.sp, bs = "re") +
                          s(elev.mean, bs = "tp", k = 8) +
                          s(elev.mean, by = bb.sp, bs = "tp", k = 8) +
                          s(fl.rich, k = 8) +
                          s(log.fl.abund, k = 8) +
                          s(bb.sp.abund, k = 5),
                        select = TRUE,
                        data = filter(spmet, year == 2011),
                        method = "REML") %>% getViz()

ggplot(filter(spmet, year == 2011), aes(elev.mean, d)) +
  geom_point() +
  geom_smooth() +
  theme_modern() +
  facet_wrap(~bb.sp)

check.gamViz(gam_d.sp_GI.2011)
concurvity(gam_d.sp_GI.2011, full = TRUE)
concurvity(gam_d.sp_GI.2011, full = FALSE)
summary(gam_d.sp_GI.2011)
print(plot(gam_d.sp_GI.2011, allTerms = TRUE), pages = 4)

gam_d.sp_GI.2012 <- gam(d ~ management +
                          s(bb.sp, bs = "re") +
                          s(elev.mean, bs = "tp", k = 8) +
                          s(elev.mean, by = bb.sp, bs = "tp", k = 8) +
                          s(fl.rich, k = 8) +
                          s(log.fl.abund, k = 8) +
                          s(log.bb.sp.abund, k = 5),
                        select = TRUE,
                        data = filter(spmet, year == 2012),
                        method = "REML") %>% getViz()

check.gamViz(gam_d.sp_GI.2012)
concurvity(gam_d.sp_GI.2012, full = TRUE)
concurvity(gam_d.sp_GI.2012, full = FALSE)
summary(gam_d.sp_GI.2012)
print(plot(gam_d.sp_GI.2012, allTerms = TRUE), pages = 4)

### I-form models
gam_d.sp_I.2010 <- gam(d ~ 
                         #management +
                         s(bb.sp, bs = "re") +
                         s(elev.mean, by = bb.sp, bs = "tp", k = 12),
                       #s(fl.rich, k = 8) +
                       #s(log.fl.abund, k = 8) +
                       #s(log.bb.sp.abund, k = 5),
                       select = TRUE,
                       data = filter(spmet, year == 2010 &
                                       bb.sp %in% c("pasc", "prat", "bss", "soro", "wurf")),
                       method = "REML") %>% getViz()

check.gamViz(gam_d.sp_I.2010)
concurvity(gam_d.sp_I.2010, full = TRUE)
concurvity(gam_d.sp_I.2010, full = FALSE)
summary(gam_d.sp_I.2010)
print(plot(gam_d.sp_I.2010, allTerms = TRUE), pages = 4)

gam_d.sp_I.2011 <- gam(d ~ 
                         #management +
                         s(bb.sp, bs = "re") +
                         s(elev.mean, by = bb.sp, bs = "tp", k = 12),
                       #s(fl.rich, k = 8) +
                       #s(log.fl.abund, k = 8) +
                       #s(log.bb.sp.abund, k = 5),
                       select = TRUE,
                       data = filter(spmet, year == 2011 &
                                       bb.sp %in% c("pasc", "prat", "bss", "soro", "wurf")),
                       method = "REML") %>% getViz()

check.gamViz(gam_d.sp_I.2011)
concurvity(gam_d.sp_I.2011, full = TRUE)
concurvity(gam_d.sp_I.2011, full = FALSE)
summary(gam_d.sp_I.2011)
print(plot(gam_d.sp_I.2011, allTerms = TRUE), pages = 4)

gam_d.sp_I.2012 <- gam(d ~ 
                         #management +
                         s(bb.sp, bs = "re") +
                         s(elev.mean, by = bb.sp, bs = "tp", k = 12),
                       #s(fl.rich, k = 8) +
                       #s(log.fl.abund, k = 8) +
                       #s(log.bb.sp.abund, k = 5),
                       select = TRUE,
                       data = filter(spmet, year == 2012 &
                                       bb.sp %in% c("pasc", "prat", "bss", "soro", "wurf")),
                       method = "REML") %>% getViz()

check.gamViz(gam_d.sp_I.2012)
concurvity(gam_d.sp_I.2012, full = TRUE)
concurvity(gam_d.sp_I.2012, full = FALSE)
summary(gam_d.sp_I.2012)
print(plot(gam_d.sp_I.2012, allTerms = TRUE), pages = 4)
```

### Does tongue-length / corolla depth matching vary with elevation?
```{r}
bb_traits <- read_csv("../processed_data/bb_traits.csv") %>%
  dplyr::select(subgenus, bb.sp, pbl.w, pbl.w.class)

fl_ktype <- read.csv("../processed_data/floral_k_type_nom_fix.csv") %>%
  dplyr::select(plant.sp, k.type) %>%
  mutate(k.type.s = str_extract(k.type, "[-+]?[0-9]*\\.?[0-9]+"),
         k.type.ss = str_extract(k.type, "[-+]?[0-9]*"))

fl_morph <- read.csv("../processed_data/flower_morphology_GitaBenadi.csv") %>%
  rename(plant.sp = Species) %>%
  group_by(plant.sp) %>%
  summarize(corolla.depth = mean(LN))

fl_tax <- read_csv("../processed_data/floral_tax.csv")

net_traits <- net %>%
  # group_by(site, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
  # summarize(freq = n()) %>%
  # ungroup() %>%
  left_join(bb_traits) %>%
  left_join(fl_ktype) %>%
  left_join(fl_morph) %>%
  separate(plant.sp, c("plant.genus", "plant.species"), sep = " ", remove = FALSE) %>%
  left_join(fl_tax) %>%
  mutate(tongue.mismatch = abs(pbl.w - corolla.depth),
         bb.sp = factor(bb.sp),
         plant.sp = factor(plant.sp),
         plant.genus = factor(plant.genus),
         subgenus = factor(subgenus),
         k.type.s = factor(k.type.s),
         plant.family = factor(plant.family))

ggplot(net_traits, aes(elev.mean, tongue.mismatch)) +
  geom_point(color = "yellow") +
  geom_smooth(method = "lm") +
  facet_wrap(~bb.sp) +
  theme_abyss()
```

### The problem with calculating network metrics on data pooled in time, particularly when you know that there is strong diversity and abundance structure in time, is that a metric like niche overlap becomes virtually uninterpretable, at least if competition is the interpretation in view. Bee 1 and Bee 2 could visit the same flower in different weeks and different flowers in the same week. In the pooled data, they might appear to have high niche overlap, but in fact they were never visiting the same flower at the same time, so there is no competition. Perhaps other metrics, like nestedness, or other interpretations -- not competition -- fare better with temporally aggregated data. So, I can consider those options. But, if I want to target an interpretation like competition through a metric like niche overlap, I have to control for time. The problem then becomes that so many of my site-dates have only a few observations, and network structure becomes uninformative (as far as I can reason) when interaction density becomes very small. You're just measuring so little of the structure that you intend to interpret. So, I've played around with the idea of setting some lower bound to drop the site-dates with too few observations, but this is inevitably arbitrary, and so many of my site-dates have very few observations that even a modest threshold culls a large number of site dates. What if I took the opposite approach, and targeted only the date(s) for each site with maximal observation density, then modeled those by elevation? The result would probably be very similar to the pooled analysis I have been doing, since pooled data becomes dominated by site-dates with high observation density. Really, all this would do is make transparent the fact that I am not controlling for time when I pool the data, I'm just masking a structure in which a few time points dominate for each site, and which time points dominate varies by site. So, this is an honest form of what I've already been doing, with no need to bootstrap.
```{r}
net.obs <- net %>%
  group_by(site, date) %>%
  summarize(total.obs = n()) %>%
  ungroup() %>%
  mutate(year = factor(year(date)))

ggplot(net.obs, aes(total.obs)) +
  geom_histogram() +
  facet_grid(site ~ year)

net.obs.top <- net.obs %>%
  group_by(site, year) %>%
  top_n(1, total.obs)

web.top <- net %>%
  semi_join(net.obs.top) %>%
  unite(site.year.date, c(site, year, date), sep = "_") %>%
  group_by(site.year.date, bb.sp, plant.sp) %>%
  summarize(freq = n()) %>%
  group_by(site.year.date) %>%
  mutate(bb.rich = length(unique(bb.sp)),
         plant.rich = length(unique(plant.sp))) %>%
  filter(bb.rich >= 2 & plant.rich >= 2) %>%
  select(higher = bb.sp, lower = plant.sp, webID = site.year.date, freq) %>%
  data.frame() %>%
  frame2webs() %>%
  map(networklevel, level = "higher", weighted = TRUE, index = c("niche overlap", "weighted NODF")) %>% # We probably care mainly about the BB level of the network
  data.frame() %>%
  t() %>%
  data.frame() %>%
  rownames_to_column() %>%
  dplyr::select(site.year.date = rowname, everything()) %>%
  separate(site.year.date, c("site", "year", "date"), sep = "_") %>%
  as_tibble() %>%
  left_join(site_data)

web.top_splevel <- net %>%
  semi_join(net.obs.top) %>%
  unite(site.year.date, c(site, year, date), sep = "_") %>%
  group_by(site.year.date, bb.sp, plant.sp) %>%
  summarize(freq = n()) %>%
  group_by(site.year.date) %>%
  mutate(bb.rich = length(unique(bb.sp)),
         plant.rich = length(unique(plant.sp))) %>%
  filter(bb.rich >= 2 & plant.rich >= 2) %>%
  select(higher = bb.sp, lower = plant.sp, webID = site.year.date, freq) %>%
  data.frame() %>%
  frame2webs() %>%
  map(specieslevel, level = "higher", index = c("d")) %>%
  map(rownames_to_column) %>%
  map(as_tibble) %>%
  bind_rows(.id = "site.year.date") %>% # Collapse list into big data frame
  rename(bb.sp = rowname) %>%
  separate(site.year.date, c("site", "year", "date"), sep = "_") %>%
  as_tibble() %>%
  left_join(site_data)

ggplot(web.top, aes(elev.mean, niche.overlap.HL)) +
  geom_point(color = "goldenrod") +
  geom_smooth() +
  theme_abyss() +
  facet_wrap(~year)

ggplot(web.top_splevel, aes(elev.mean, d)) +
  geom_point(color = "goldenrod") +
  geom_smooth() +
  theme_abyss() +
  facet_wrap(~year)
```


### Play with null models
```{r}
webs <- net %>%
  group_by(site, year, bb.sp, plant.sp) %>% # group by site, date, bee species, plant species
  summarize(freq = n()) %>% # calculate interaction frequency per species pair within each site
  ungroup() %>%
  unite(webID, c(site, year), sep = "_") %>%
  dplyr::select(higher = bb.sp, lower = plant.sp, webID, freq) %>% # rename columns to match bipartite's expectations
  data.frame() %>% # convert to data frame
  frame2webs()

null.t.test(webs[[2]], level = "higher", index = "niche overlap")

n1 <- nullmodel(webs[[1]], N = 1000, method = 1)
met1 <- map(n1, specieslevel, level = "higher", index = "resource range") %>%
  bind_rows(.id = "id") %>%
  
  ```