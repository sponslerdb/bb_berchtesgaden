---
title: "Bumble bees and wildflowers along a mountain slope exhibit nonlinear patterns of abundance and turnover, punctuated by the tree line ecotone"

author: Douglas B. Sponsler^1^*, Fabrice Requier^2^, Katharina Kallnik^1^, Alice Claßen^1^, A. Fabienne Maihoff^1^, Ingolf Steffan-Dewenter^1^

# date: 2021-02-18

output:
  bookdown::word_document2:
    fig_caption: yes
    fig_width: 7
    fig_height: 7
    df_print: kable
    keep_md: yes
    number_sections: no
    toc: no
    reference_docx: Word_template.docx
    
bibliography: 
  - references.bib
#  - knitcitations.bib
csl: ecology.csl
# csl: https://www.zotero.org/styles/ecology

---

^1^ Department of Animal Ecology and Tropical Biology, Biocenter, University of Würzburg, Würzburg, Germany

^2^ Université Paris-Saclay, CNRS, IRD, UMR Évolution, Génomes, Comportement et Écologie, 91198, Gif-sur-Yvette, France

`*` Corresponding author: douglas.sponsler@uni-wuerzburg.de


# ABSTRACT

Write your abstract here.


*Keywords*: climate, elevation gradient, pollinator, *Bombus*, alpine









# INTRODUCTION

The rate of modern climate change is forcing species to adapt ecologically within the constraints of gene pools shaped by cooler climates [@Visser2008-ct]. As temperatures rise, species are expected to track suitable conditions to higher latitudes and higher elevations, and such effects are already measurable across a broad range of taxa [@Lenoir2015-tj]. Species range shifts, however, are not merely abstract, Cartesian translations along axes of latitude, elevation, and temperature; they are movements *in situ*, involving the advent, extinction, and reorganization of biological interactions [@Tylianakis2008-fi; @Blois2013-dz] as well as the concomitant variation of abiotic conditions other than temperature [@Korner1995-mh; @Hodkinson2005-af]. Moreover, some of these factors, unlike temperature, cannot be assumed to vary linearly with latitude or elevation but instead exhibit marked nonlinearities and thresholds. A salient example is the tree line ecotone that marks the elevational transition between forest and grassland. At this relatively discrete band within a continuous elevational temperature gradient, plant communities exhibit rapid species turnover [@Descombes2017-id] and declining diversity [@Becker2007-cg] as abiotic conditions become more extreme away from the shelter of tree canopy [@Slatyer1992-ir]. Thus, species range shifts --- and, perhaps more importantly, the ecological consequences thereof --- are contingent upon patterns of species distributions and biotic interactions, and community ecology is as much prior as posterior to the ecology of global change. 

Bumble bees (Hymenoptera: *Bombus* spp.) are both a classic model organism in community ecology and an emergent conservation priority due to recent declines [@Goulson2008-oe] and manifest vulnerability to climate change [@Soroye2020-zb]. The genus *Bombus* is believed to have arisen in the mountains of Asia during the cooling climate of the Eocene-Oligocene boundary (~34 mya), and its subsequent spread and diversification along mountain corridors and into lowland habitats was likely driven by alternating range expansion during cooling periods and retreat to higher elevations during warming periods [@Hines2008-er; @Stewart2010-jn; @Martinet2018-rb]. Today, bumble bees remain characteristically cold-adapted species [@Heinrich1994-cv], with peak abundance and diversity associated with mountain ranges and northern latitudes, where they are often the dominant guild of flower-visiting insects and the principal pollinators of entomophilous flora [@Williams1998-xi; @Goulson2010-lo]. Climate-induced range shifts are expected to be especially pronounced in the mountain habitats that host the core of bumble bee biodiversity, together with their floral mutualists, and upward movement of bumble bee species has already been documented [@Ploquin2013-da; @Kerr2015-dg; @Pyke2016-ll; @Fourcade2019-ct; @Soroye2020-zb; @Marshall2020-pc]. Populations of bumble bee species already restricted to the highest elevation zones can be expected to suffer geometric habitat shrinkage and increasing genetic isolation as they move toward mountain peaks; the upward advance of lowland species may induce novel competitive pressures and disrupt the pollination networks in which bumble bees participate [@Brosi2013-mu; @Ishii2013-sm]; and concurrent elevational [@Lenoir2008-yn] and/or phenological [@Pyke2016-ll] shifts in the floral community could generate mismatches between bumble bee species and their historic floral hosts. 

The outcomes of range shifts in mountain bumble bee and floral communities have obvious significance from a conservation perspective, but they also provide a model system in which to study the general phenomenon of community assembly under climate change, a phenomenon of both historical interest and prospective urgency. The goal of anticipating, interpreting, and potentially mitigating the effects of ranges shifts, however, is predicated upon understanding the current elevational and temporal patterns of bumble bee communities and those of their floral hosts. In particular, it is important to detect nonlinearities in species-, community-, and interaction-level responses to elevational climate gradients, such as the floristic threshold effect that has been documented at the tree line ecotone [@Descombes2017-id]. 

In the present study, we approach this topic at three levels of biological organization. First, we ask how mountain bumble bee species are distributed in elevation and intra-annual time with respect to one another and with respect to their floral hosts. Then, we move from the level of individual species to that of community composition and explore the structure of species turnover (i.e. β-diversity) of bumble bees and flora in response to elevation. Finally, propagate the question to the level of species interactions, investigating the turnover of bee-flower interaction partners through elevation. In each of these analyses, we employ flexible modeling techniques capable of capturing nonlinear responses to elevation and time with the goal of identifying where, when, and for which species the impacts of climate-induced range shifts should be expected to be most acute.  



# METHODS

## Field system

The study was conducted in Berchtesgaden National Park (47.55°N, 12.92°E), located in the Northern Limestone Alps of southeast Germany. The landscape is composed of mountain pastures mainly surrounded by coniferous forests. We selected 25 study sites (60 x 60 m) on mountain pastures at elevations ranging from 641-2032 m above sea level. Fourteen of these pastures are extensively grazed by cattle or sheep, three are mowed for hay production, and eight have been abandoned throughout the last century and are no longer subject to any human management. The sites sampled in this study are the same as those used in a series of previous studies [@Hoiss2012-sn; @Hoiss2012-kr; @Hoiss2015-gv], but the data reported here are otherwise independent.

Sampling consisted of repeated visits to each study site at approximately weekly intervals. Samples were only collected on dry days when the air temperature was at least 6°C. During each visit, bumble bees and their activity were recorded during a 50-minute transect walk. Bumblebees observed on or a given flower were counted as floral visitors. Bumble bee queens were identified to species level in the field, while workers and males were stored in individually labeled tubes in the freezer for later identification in the laboratory after [@Amiet1996-jz]. Floral visitation by male bumble bees was recorded during visitation sampling, but we chose to analyze only visitation by queens and workers. 

In conjunction with visitation observations, we estimated the flower cover of each herbaceous or shrubby plant species within each 60 x 60 m study plot to the nearest 0.1 m^2^. Species identification followed @Lauber2007-ha and @Oberndorfer2001-ro.


## Data analysis

Bumble bee abundance was quantified as the total number of recorded floral visits per bumble bee species per site-date (excluding males). To represent floral abundance from the perspective of each bumble bee species in our network, we first scored each plant genus in our data set as either visited or not visited by each bumble bee species. We then calculated the total floral resource availability for each bumble bee species for each site-date by summing the observed flower cover of visited genera. We opted to work at the genus level based on the reasoning that if a bumble bee species visits one member of a given genus, other members of the same genus should also be considered potential floral hosts. This approach was intended to dampen the effects of false non-detection on the estimation of floral resource abundance, particularly for rare bumble bee species whose diet breadth would tend to be underestimated simply due to sparsity of observations [@Williams2005-yp].

We analyzed bumble bee abundance and floral abundance using separate hierarchical generalized additive models (HGAM) [@Wood2011-pg; @Pedersen2019-kf], with abundance modeled as a response to the interaction (tensor product) of elevation and day-of-year. Bumble bee species was included in each model as both a smoothing factor and a random intercept term, resulting in centered 2-dimensional abundance smooths for each bumble bee species. We also included year and site as random intercept effects to account for annual differences in abundance patterns and repeated measurements within sites. 

To complement GAM-based analysis of abundance patterns through elevation and time, we also analyzed the site occupancy pattern of each bumble bee species using cluster and ordination analyses. We first calculated the mean abundance across dates of each bumble bee species for each site in our study system. Then, we applied a Hellinger transformation to the matrix of mean abundances and calculated the Horn-Morisita distance metric for each pair of species. The resulting distance matrix was then used to generate a principal coordinates analysis (PCoA) and a cluster dendrogram, each depicting the (dis)similarity between species pairs in their distribution across sites. The primary purpose of this analysis was to organize the presentation of the GAM results by grouping species with similar site affinities, but because our sites were arranged along an elevation gradient, our cluster and ordination analyses also highlight differences between species in elevational distribution.

For our analyses of species and interaction β-diversity, we opted to focus on elevational patterns by pooling samples within sites, since the sparse observations of individual sampling events tended to inflate β-diversity and its variance to a degree that impaired model fitting and interpretation. After pooling, we calculated the species and interaction β-diversity between all pairs of sites and partitioned total interaction β-diversity into its components of species turnover (further partitioned into bumble bee and plant turnover, respectively) and interaction rewiring [@Novotny2009-my; @Poisot2012-fk]. 

As an initial analysis, we plotted the relationship between site-wise elevation difference and each metric of β-diversity and verified the significance of the relationship using logistic matrix regression [@Goslee2009-ln], with the difference in number of sampling dates between sites as a covariate to control for potential confounding effects of sampling frequency.

To investigate the relationship between elevation and β-diversity more deeply, we performed a second analysis using generalized dissimilarity modeling (GDM) [@Ferrier2007-ji] to analyze the elevational variation in the β-diversity of bumble bees, flora, and (unpartitioned) interactions. In addition to estimating the relationship between differences in the response variable and differences in the predictor variables, GDM captures the slope of this relationship over the range of each predictor variable, revealing potential variation in the amount of change in the response induced by a given change in a predictor. This enabled us to ask whether the relationship between β-diversity and elevation exhibits thresholds or other nonlinearities. To control for potential confounding effects of geographic proximity and number of sampling dates per site, these terms were added as covariates in the GDM.

All analyses were conducted in `R` [@R_Core_Team_2020]. GAM analyses were performed with packages `mgcv` [@Wood_2017] and `mgcViz` [@Fasiolo2018aa], and the models we fit correspond to the “type I” model form described in @Pedersen2019-kf. PCoA and species β-diversity calculation were performed with the package `vegan` [@Oksanen_2019]. Calculation and partitioning of interaction β-diversity were performed with the package `bipartite` [@Dormann_2009; @Dormann_2008; @Dormann_2011]. Matrix regression was performed with the package `ecodist` [@Goslee2007-aa]. GDM analysis of β-diversity was performed with the package `gdm` [@Fitzpatrick_2020]. Data handling and visualization was performed with the `tidyverse` suite of R packages [@Fitzpatrick_2020]. Annotated R code is available in the Supplementary Material. 


# RESULTS

## Summary of visitation and floral survey data

We recorded a total of 12,918 bumble- bee-flower interactions (excluding males) over the three years of our study. The metaweb across all sites and dates consisted of 16 bumble bee species (with *Bombus terrestric/lucorum* and *Bombus psithyrus* species groups lumped into species groups), 163 plant species (110 genera, 37 families), and 736 unique bumble- bee-plant interaction pairs. Five species --- *B. pascuorum*, *B. pratorum*, *B. soroensis*, *B. terrestris-lucorum*, and *B. wurflenii* --- accounted for the bulk of overall bumble bee abundance and were present over the entire elevation range of our sites **(Figure 1)**. *B. hortorum*, *B. jonellus*, and *B. psithyrus* (species-group) were widespread but at lower abundance. *B. monticola* occurred at moderate abundance but was rare or absent at sites below 1000 m. The remaining species were relatively rare in our study system. *B. mendax*, *B. mucidus*, and *B. pyrenaeus* occurred primarily above 1500 m, consistent with their known affinity for high elevation habitats [@Rasmont2010-ci]. *B. gerstaeckeri*, *B. hypnorum*, *B. lapidarius* were sparsely distributed across mid-elevation sites. *B. humilis* was recorded only three times and only in one year, so we omitted it from all analyses.

![(\#fig:fig1)Mean abundances and elevational distributions of bumble bee species. Sampled points along the elevation gradient are depicted as circles, the diameter of which is scaled to the mean abundance (pooled across study years) of the given bumble bee species at that elevation. Abundances of zero are not plotted. Lines span between the highest and lowest occurrence of each species, highlighting its observed elevation range.](../output/bb_range.png)

Floral surveying yielded a total of 352 plant species, representing 191 genera and 52 families. Of these, 155 species from 103 genera and 35 families were observed to be visited by bumble bees. Eight species --- *Rubus idaeus*, *Rosa canina*, *Juniperus communis*, *Larix decidua*, *Salix* sp., *Caltha palustris*, *Pulmonaria officinalis*, and *Rheum barbarum* --- were recorded in visitation data but not recorded during floral surveying. Each accounted for no more than 4 visits in total over the three years of our study, and they were omitted from the analysis of floral abundance.

## Abundance over elevation and time

Cluster **(Figure 2A)** and ordination **(Figure 2B)** analyses highlighted the strong differences in site affinity between the high-elevation species *B. mucidus*, *B. monticola*, *B. mendax*, and *B. pyrenaeus* and the rest of the bumble bee community. PCoA captured 59% of total variation in site affinity. 

Elevation-time GAM smooths **(Figure 2C)** of bumble bee abundance differed significantly from a flat surface (p < 0.00001) for all species except *B. hypnorum* (p = 0.8). Overall, the model explained 57.6% of total deviance (adjusted R^2^ = 0.346). The distinction between *B. mucidus*, *B. monticola*, *B. mendax*, and *B. pyrenaeus* and the rest of the community was again clear; all 4 high-elevation species exhibited peak abundance at the upper extreme of our elevation gradient, and only *B. monticola* was found at high abundance below 1500 m. *B. monticola* also peaked in abundance earlier in the season than the other high-elevation species. 

*B. gerstaeckeri* and *B. pratorum* exhibited well-defined mid-elevation abundance peaks, with that of *B. gerstaeckeri* occurring very late in the season. The apparent mid-elevation abundance peak of *B. hypnorum* was not significant and should not be interpreted strongly. 

*B. pascuorum* and *B. psithyrus* had well-defined low-elevation peaks, with that of *B. psithyrus* occurring early in the season, during the colony-founding --- or, in the case of *B. psithryus*, colony-usurping --- phase of the bumble bee life cycle. 

The remaining species exhibited relatively diffuse patterns of abundance. Of these, *B. soroensis* was exceptional in that it exhibited peak abundance at the highest sites, but it was also found at relatively high abundance across the whole elevation gradient, peaking at mid-season. *B. hortorum* and *B. lapidarius* peaked at low elevation. *B. jonellus* exhibited little elevational variation in abundance until declining above 1500 m, but its abundance was strongly patterned in time, peaking early in the season. *B. terrestris/lucorum* was generally abundant across the whole elevation gradient and throughout the whole season, but it showed a weakly trimodal pattern in elevation, with abundance peaks at low, middle, and high sites. Its strongest abundance peak occurred at low elevation, early in the season during colony-founding. *B. wurflenii* was broadly abundant above 1000 m but was rarer at lower sites.

Elevation-time GAM smooths of floral abundance differed significantly from a flat surface for all species (p << 0.00001) and explained 64.8% of total deviance (adjusted R^2^ = 0.425). Floral resource availability exhibited similar patterns for most bumble bee species, with peak floral resources occurring below 1500 m and during mid-season **(Figure 2)**. There was a noticeable temporal lag in floral resource availability associated with increasing elevation, consistent with expected temperature constraints on plant phenology. The only species to exhibit marked anomalies in floral resource availability was *B. gerstaeckeri*, which specializes on the genus *Aconitum*[@Ponchau2006-af] but was also observed (presumably collecting nectar) on *Gentiana*, *Salvia*, *Rhinanthus*, *Cirsium*, and *Clinopodium.* Its floral resources were concentrated in two discrete peaks, the first at low elevation and mid-season, and the second at high elevation and late-season. 


![(\#fig:fig2)Abundance of bumble bee species and their floral resources through elevation and time. For each bumble bee species, surface plots are arranged in two columns labeled BB and FL. Surface plots in the BB columns show the abundance of the bumble bee species, while plots in the FL columns show the abundance of floral resources for the corresponding bumble bee species. The color ramp of each smooth has been scaled to equal color range, thereby obscuring the absolute magnitude of effect. To avoid cluttering the figure, we have opted not to include a unique color ramp key for each plot, since the focus of our inference is on relative rather than absolute patterns. Surface plots with original color ramp scales are available in Appendix X.](../output/gam_BB_FL4.png)

## Species and interaction β-diversity

Floral β-diversity among our study sites was very high overall and responded steeply to elevation difference between sites **(Figure 3A)**. Even sites at similar elevation exhibited ~40% species turnover, and the most widely separated sites (elevation difference > 1250 m) differed by more than 85%. 

Bumble bee β-diversity was, in comparison to floral β-diversity, both lower overall and less responsive to elevational difference between sites **(Figure 3A)**. Sites at similar elevation exhibited ~25%
species turnover, and species turnover between the most widely separated sites remained less than 50%. 

Total interaction β-diversity was >75% between sites at similar elevation and approached perfect dissimilarity in the most widely separated sites **(Figure 3B)**. Partitioning revealed that total interaction β-diversity was driven primarily by species turnover, which accounted for ~50% of total β-diversity for sites at similar elevation and ~90% of total β-diversity for the most widely separated sites. Species turnover consisted mainly of turnover in the floral community **(Figure 3C)**, though the joint turnover of plants and bumble bees accounted for ~25% of total species turnover between the most widely separated sites. Interaction rewiring was most significant between sites at similar elevation, where it accounted for ~30% of total β-diversity, but its share of total β-diversity declined to ~10% in the most widely separated sites **(Figure 3B)**.

Logistic matrix regression confirmed the significance (p < 0.01) of the relationship between each β-diversity metric and elevation **(Table 1)**.

![(\#fig:fig3)Species-level (A) and interaction-level (B, C) plotted against elevational difference. Each point represents the β-diversity (or β-diversity partition) between a pair of sites, and lines represent the overall relationship between β-diversity and elevation difference using binomial regression smooths. Standard errors are not plotted because they would be misleading due to the non-independence inherent to distance matrix regression, but all regressions were significant (p < 0.01). Interaction β-diversity (B, C) is partitioned using Poisot et al.’s (2012) notation: WN = unpartitioned β-diversity, ST = β-diversity due to species turnover, OS =  β-diversity due to interaction rewiring, ST.h =  β-diversity due to species turnover in the higher trophic level (bumble bees), ST.l = β-diversity due to species turnover in the lower trophic level (plants), and ST.lh =  β-diversity due to joint species turnover in higher and lower trophic levels.](../output/beta_stack.png)

GDM analysis reproduced the finding that overall β-diversity was dominated by floral turnover, but it also revealed a strongly nonlinear response of floral and interaction β-diversity to elevation **(Figure 4)**. Turnover was high from ~600-1000 m, representing the transition from valley floor to lower slopes, and then leveled off from ~1000-1500 m. At around 1500 m, turnover accelerated sharply and remained steep for the remainder of the elevation gradient. The partial effects of geographic proximity and sampling intensity were negligible in all cases.

![(\#fig:fig4)GDM splines of species-level and interaction-level β-diversity in response to geographic, elevational, and temporal distance between site-dates. The maximum height of each spline represents a given variable’s partial effect on β-diversity (its effect when covariates are held constant), and the shape of each spline represents to the rate of species turnover as it varies along the corresponding gradient (i.e. geographic distance, elevation, or day of year), with steeper parts of the curve indicating regions of the gradient over which species/interaction turnover is more rapid. Within-site comparisons were omitted from GDM analysis to avoid pseudoreplication. Shaded bands depict the uncertainty of each partial effect curve (+/- one standard deviation) based on bootstrapping. Note the differences in scale between the y-axes of each variable.](../output/gdm_stack.png)

# DISCUSSION



The floristic importance of the tree line ecotone has been noted in previous work [@Pellissier2010-kg; @Descombes2017-id]. By studying the floral community through the lens of foraging bumble bees, we provide a functional extension of this floristic pattern, showing that the tree line ecotone both a boundary above which floral resource availability sharply declines and a critical distributional interface between high-elevation bumble bees and their habitat-generalist counterparts. Similar patterns with respect to the tree line ecotone have been described for bumble bees and their floral hosts in the Mediterranean system of Mount Olympus [@Minachilis2020-xd], suggesting that these findings are a general property of mountain ecosystems. 

The fact that

As such, the tree line ecotone deserves special consideration in conservation management, and targeted studies of plant-pollinator relationships at the tree line ecotone would be warranted.

A salient question that could not be answered by our sampling approach is whether mountain bumble bees adaptively forage up- or downslope, as suggested by @Lundberg1980-qk. Bumble bees have large foraging ranges and have been shown both to track resources through space and time [@Devoto2014-ay] and to cross forest matrix to reach patches of foraging habitat [@Mola2020-yc]. Assuming an average slope of 20° and a foraging range of 1 km, a bumble bee could travel up- or downslope by more than 340 m, thus spanning an elevation belt nearly 700 m wide. Such 3-dimensional foraging would enable bumble bees both to exploit the elevational turnover of floral species and to track preferred species through their elevationally staggered phenology, the latter constituting a sort of physiological time travel [@Van_Straalen1983-df. 

he rate of β-diversity --- whether that of bumble bees, flora, or interactions --- accelerates sharply near the tree line ecotone (~1500 m) and continues accelerating for the remainder of the elevation gradient (Figure 8). Thus, a linear elevation gradient generates a nonlinear biological response, punctuated by the tree line ecotone, above which further increase in elevation becomes an increasingly stringent filter of community assembly. Such patterns have been described previously for plants (Pellissier et al. 2010, Descombes et al. 2017); our study suggests they extend also to bumble bees, and that they are propagated from the level of community composition to that of species interactions between bumble bees and their floral hosts.


# CONCLUSIONS

The central finding of our study, evident in each of our analyses, is that a linear elevation gradient can generate nonlinear biological responses in terms of both abundance and β-diversity, and the latter at both at the level of community composition and at the level of species interactions. While elevational climate variation cannot be completely disentangled from other variables that covary with elevation [@Hodkinson2005-af], nonlinear responses of plants and pollinators along an elevational climate gradient suggest that we might expect similar nonlinearities in response to temporal climate change. Indeed, such nonlinear responses to climate have already been reported, including a recent study that documented nonlinear decline of insect abundance over a 37-year temporal climate gradient [@Lister2018-jp].

Independent of their utility as space-for-time climate proxies, mountains per se are unique ecological theaters in which extrinsic temporal climate change interacts with an intrinsic elevational climate gradation to alter the spatial and temporal distributions of species and interactions [@Telwala2013-sc; @CaraDonna2014-nq; @Miller-Struttmann2014-wz; @Rafferty2020-aj]. The nonlinearity of abundance observed in our study warns that temporal climate change on mountain slopes should be expected to involve the redistribution of discrete peaks of species abundance.


It should also be expected that climate change on mountain slopes will be punctuated by the tree line ecotone, which constitutes, at least with respect to bumble bees and wildflowers, a region of unique resource abundance and a transition to high and accelerating turnover of species and interactions [@@Descombes2017-id]. The observation that the upslope advancement of alpine tree lines often lags behind climate warming [@Dullinger2004-qu] suggests a danger that high elevation florivores could move upslope faster than their floral resources, or perhaps be forced to cope with increased competitive pressure from other florivore species formerly restricted to lower elevations [@Hoiss2012-sn]. 

though as a matter of applied conservation it is perhaps doubtful that local interventions could forestall for long the consequences of a process occurring at a global scale. 

# ACKNOWLEDGEMENTS

We thank D.J. McNeil for helpful conversations during data analysis.


# REFERENCES

<!-- ```{r write_citations, cache=FALSE, include=FALSE} -->
<!-- write.bibtex(file = "knitcitations.bib") -->
<!-- ``` -->

<div id ="refs"></div>
