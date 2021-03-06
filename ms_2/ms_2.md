---
title: "Bumble bees exhibit strong but flexible foraging biases through temporal and elevational gradients of floristic turnover"

author: Douglas B. Sponsler^1^*, Fabrice Requier^2^, Katharina Kallnik^1^, Alice Claßen^1^, A. Fabienne Maihoff^1^, Ingolf Steffan-Dewenter^1^

# date: 2021-02-25

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



*Keywords*: optimal foraging, trait matching, pollinator, *Bombus*, alpine









# INTRODUCTION

<!-- Open with a discussion of the problem of foraging efficiently in a dynamic foraging environment, the balance of bias and flexibility, etc. -->

Bumble bees are exceptionally versatile foragers whose ability to solve the puzzles posed by complex flower morphologies (Plowright and Laverty 1984), or indeed by human investigators (Alem et al. 2016; Loukola et al. 2017), renders almost any flower potentially accessible. Even when a gross mismatch between tongue length and corolla depth makes the nectar of a given flower inaccessible by ordinary means, some bumble bee species will resort readily to biting through the base of the flower to “rob” nectar, incidentally making the nectar available to other prospective foragers as well. Flower handling skills are, moreover, socially transmitted, amplifying the plasticity of individual behavior (Leadbeater and Chittka 2008; Loukola et al. 2017). Nevertheless, the flexibility of bumble bee foraging is tuned to floral reward, which depends on both the intrinsic reward of a given flower and the searching/learning/handling costs associated with exploiting it. At an evolutionary level, the innate biases of bumble bees for certain floral traits likely reflect historic associations between floral traits and floral rewards (Raine and Chittka 2007). These innate biases, however, can be modulated or even nullified by patterns of reward learned during foraging (Kugler 1943; Maharaj et al. 2019), which depend not only on intrinsic floral properties but also ecological modulators, such as the spatial distribution of floral patches or competition with other bumble bees (Ishii 2013; Plowright and Laverty 1984; Sowig 1989).

The topic of floral bias in bumble bees has been approached via both laboratory experimentation (e.g. Kugler 1943; Raine and Chittka 2007) and observational field study (e.g. Heinrich 1976; Dave Goulson, Lye, and Darvill 2008). The former approach enables stronger causal inference with respect to the mechanisms of floral choice, but it suffers from all the usual limitations of biological realism associated with laboratory conditions. Observational field studies, however, tend to be anecdotal in scope, typically describing bumble bee foraging patterns in particular localities, over short time frames, and with limited floristic data to enable a distinction between true foraging bias and the neutral effects of floral abundance. The steep elevation gradients of mountain slopes provide an opportunity for more powerful observational study of bumble foraging because they generate broad gradients of environmental variation while controlling for confounding effects of geographic variation, such as genotypic differentiation within bumble bee species. When bumble bee visitation in such environments is sampled with sufficient spatial and temporal replication and complemented with quantitative floral surveying, true bias can be distinguished from the neutral effects of abundance, and variation in bias can be explored along elevational and temporal dimensions of floristic turnover.

In the present study, we focus on the question of bumble bee foraging bias with respect to floral morphology, and the lability of morphological biases across gradients of floristic variation. At the outset, we emphasize the use of the term “bias” in favor of the term “preference,” since the latter is laden with more assumptions and its inference under field conditions is problematic (Williams 2005). By “bias,” we mean simply an empirical disproportionality between the relative abundance of a given floral taxon or morphological class and its proportional visitation. The exact mechanism of this disproportionality --- whether an arbitrary “preference” on the part of the foraging bee, an unmeasured decoupling of relative abundance and rate of reward, an exogenous effect of competition, etc. --- is beyond the scope of our discussion and probably best studied under laboratory conditions where strict factorial control is possible. We focus instead on the empirical phenomenon of bias, its variation, and its morphological correlations.

The meaning of “morphology,” however, is hardly more clear-cut than that of bias. Bumble bees and flowers vary in an array of interrelated morphological parameters, both within and between species, and the selection of relevant and measurable morphological metrics is nontrivial. For bumble bees, morphological variation consists principally (though not exclusively) of variation in tongue length. Tongue length in bumble bees corresponds functionally to corolla depth in flowers, and there is a well-established (though imperfect) association between tongue length and floral choice (e.g. Inouye 1980). Importantly, though, accessing floral reward involves more than just a match between tongue length and corolla depth. Such matching is entirely bypassed by nectar robbers, but even in the case of “legitimate” foraging, access to reward involves learning to manipulate floral structures that, in some floral morphotypes, are challenging irrespective of tongue length (Laverty 1994). For example, the characteristic Papilionaceous floral morphology of legumes requires a bee to separate tightly appressed keel petals to access anthers and nectaries (Westerkamp 1997), and the narrowly campanulate flowers of many Ericaceous flora pose a significant barrier to large-bodied visitors like bumble bees (Wallace 1977). Rather than considering various floral morphological traits independently, we rely on the synthesis of Kugler (1970), who developed a typology consisting of ten primary morphotypes and nested subdivision of each. Kugler’s typology was informed by laboratory studies of bumble bee flower choice (Kugler 1943), and there is a precedent for analyzing bumble bee visitation in terms of Kugler’s morphotypes (Neumayer and Paulus 1999; Schneller et al. 2014).


<!-- Formulate this section more as a list of motivating questions; perhaps end with something that ties back into the more general introduction. -->

How do floral morphotypes compare in terms of total bumble bee visitation?

How does morphological bias vary across bumble bee species and tongue length classes?

Species' morphological biases give rise to patterns of diet (dis)similarity 

How are the floral morphotypes distributed through elevation and time in our mountain ecosystem?

For ubiquitously distributed bumble bee species, how does morphological bias vary through elevation and time?

First, we compare the overall rate of bumble bee visitation across floral morphotypes. Then, we parse visitation rates by bumble bee species and explore the functional relationship between floral morphotype and bumble bee tongue length. These biases, together with the background effects of relative abundance, generate emergent patterns of diet (dis)similarity, which we present in the form of ordination analysis. Lastly, we analyze the lability of the morphological biases of individual bumble bee species through elevation and time, and compare these patterns of bias across bumble bee species and within floral morphotypes. The results of these analyses demonstrate both the strength of bumble bee floral bias, arising in part from morphological constraints, and the flexibility of bumble bee floral bias that enables them to forage adaptively in the dynamic floral communities of mountain landscapes.
  



# METHODS

## Field system

<!-- I guess I'll need to tweak the wording of this section so that it is not a verbatim copy of the other manuscript. But it's okay if it recycles a lot of language because it is, after all, the same thing being described. -->

<!-- Describe general floristic and temperature patterns. Describe sampling in more detail. -->

The study was conducted in Berchtesgaden National Park (47.55°N, 12.92°E), located in the Northern Limestone Alps of southeast Germany. The landscape is composed of mountain pastures mainly surrounded by coniferous forests. We selected 25 study sites (60 x 60 m) on mountain pastures at elevations ranging from 641-2032 m above sea level. Fourteen of these pastures are extensively grazed by cattle or sheep, three are mowed for hay production, and eight have been abandoned throughout the last century and are no longer subject to any human management. The sites sampled in this study are the same as those used in a series of previous studies [@Hoiss2012-sn; @Hoiss2012-kr; @Hoiss2015-gv], but the data reported here are otherwise independent.

Sampling consisted of repeated visits (hereafter, “visitation samples”) to each study site at approximately weekly intervals. Samples were only collected on dry days when the air temperature was at least 6°C. During each visit, bumble bees and their activity were recorded during a 50-minute transect walk. Bumblebees observed on or a given flower were counted as floral visitors. Bumble bee queens were identified to species level in the field, whereas workers and males were stored in individually labeled tubes in the freezer for later identification in the laboratory after Amiet (1996). While floral visitation by male bumble bees was recorded during visitation sampling, we chose to analyze only visitation by queens and workers. During observation, it was noted whether each bee was visibly carrying pollen at the time of its visit, but this provides at best an imperfect distinction between pollen and nectar foragers, and we opted not to differentiate between bees with and without pollen in our analyses. In conjunction with visitation observations, we estimated the flower cover of each herbaceous or shrubby plant species within each 60 x 60 m study plot to the nearest 0.1 m^2^. Species identification followed @Lauber2007-ha and @Oberndorfer2001-ro.


## Data analysis

### Data processing

Bumble bee tongue-length data were obtained from the literature values compiled and summarized by Arbetman et al. (2017). Values for *B. gerstaeckeri* and *B. mendax*, which were not included in Arbetman et al. (2017), were obtained from Obeso (1992) and Durieux (2000), respectively. We then binned bumble bee species into three tongue-length classes: short (< 7.5 mm), medium (7.5—10.0 mm), and long (>= 10.0 mm). Short-tongued species included *B. jonellus*, *B. lapidarius*, *B. pratorum*, *B. psithyrus*-group, *B. pyrenaeus*, and *B. terrestris/lucorum*. Medium-tongued species included *B. humilis*, *B. hypnorum*, *B. monticola*, *B. mucidus*, *B. pascuorum*, and *B. wurflenii*. Long-tongued species included *B. gerstaeckeri*, *B. hortorum*, and *B. mendax*.

Floral morphology was summarized using floral morphotypes of Kugler (1970), accessed via the BIOFLOR database (Klotz, Kühn, and Durka 2002). Morphotypes were simplified to the primary classes numbered 0-10: 0 = not applicable (includes wind-pollinated flora such as grasses), 1 = disc- and bowl-shaped flowers, 2 = funnel flowers, 3 = bell-shaped flowers, 4 = stalk-disc flowers, 5 = lip flowers, 6 = flag flowers, 7 = flower heads, 8 = spike flowers, 9 = brush flowers, 10 = trap flowers **(Figure 1)**. In a few cases, a plant species in our data set was not found in the BIOFLOR database. In most cases, morphotypes are constant at the genus level, and missing morphotypes were added based on the morphotypes found for congeners. The genus *Gentiana* consisted of both group 2 and group 4 representatives, so the missing species *Gentiana aspera* and *Gentiana ciliata* were classified based on visual comparison with their classified congeners. In the case of the genus *Phyteuma*, we followed Neumayer and Paulus (1999) and classified it as a group 7 flower head.

Prior to analysis, we aligned each visitation sample with a corresponding floral survey sample. In most cases, floral surveying was conducted jointly with bumble bee observation for each site and date. Occasionally, though, floral surveying could not be completed on the same day as bumble bee observation, resulting in a small temporal offset between visitation samples and floral survey samples. In these cases, visitation samples were paired with the nearest available floral survey data within a maximum distance of 7 days. Visitation samples lacking corresponding floral survey data within 7 days were omitted from further analysis. In a few cases, two consecutive visitation samples at a given site shared the same nearest floral survey sample; in these cases, the visitation samples were first pooled, then aligned to their shared nearest floral survey sample. In rare cases, a given floral species was recorded during visitation sampling but missed during floral surveying. These species were retroactively assigned the minimum floral area value of 0.01 m^2.

Visitation data were extremely right-skewed because the majority of available flora received zero visits during a given sampling period. To facilitate model fitting, we created a binary version of visitation rate in which visited flora were assigned to 1 and unvisited flora to 0. We then modeled binary visitation using generalized models with binomial or Bernoulli distribution families. The conditional effects inferred from these models should, therefore, be interpreted as effects on visitation probability rather than visitation frequency.

<!-- Should I be including year as a varying effect in these models? -->

<!-- State more explicityl how the models should be interpreted withr expect to bias -->

<!-- Report the effects of abundance that I am conditioning on -->

<!-- Perhaps show what happens if we remove the abundance interaction terms -->

<!-- Mor edescriptive fgure captions -->

<!-- Add something about model diagnostics -->

<!-- Consider adding year as a varying effect -->

<!-- Figure 3 order by tongue length class -->

<!-- Describe the differences in models in the Figure2 and 3 captions -->

<!-- Explain filtering of ordination plot -->

<!-- Year should really be included in Model 5, or just use 2012. -->

<!-- Open each results paragraph with a summary sentence. -->

### Model 1: Comparing overall visitation probability across floral morphotypes

We first compared overall visitation probability across floral morphotypes. In this model, visitation for each sample was pooled across bumble bee species and within floral morphotypes prior to binarization. Visitation probability was then modeled using multilevel Bayesian regression with a Bernoulli distribution family (logit link function). *Bumble bee abundance* and *percent flower cover* (per morphotype) were included as a global interaction term, thus conditioning visitation probability on abundance and enabling the inference of bias. *Floral morphotype* was specified as a varying effect term with both a varying intercept and a varying slope with respect to *percent flower cover*, and a varying intercept term was also included for *site*. See Appendix X for the details of model specification and evaluation.

### Model 2: Visitation probability as a function of floral morphotype and bumble bee tongue length

To explore the functional underpinnings of bumble foraging bias, we analyzed visitation probability as a function of the interaction between floral morphology and bumble bee tongue length class. For this analysis, visitation was pooled within floral morphotype but not across bumble bee species. We then specified a Bayesian model of the same form as Model 1, but this time included a varying slope and intercept term for floral morphotype : bumble bee tongue length class and floral morphotype : bumble bee species. In this model, the relationship between floral morphotype and bumble bee tongue length class was the focus of inference, and the bumble bee species effect was added to avoid pseudoreplication arising from the grouping of observation within species. See Appendix X for the details of model specification and evaluation.

### Model 3: Visitation probability as a function of floral morphotype and bumble bee species 

To focus on species-level morphotype biases, we removed the tongue length effect from Model 2 and re-fit the model considering only bumble bee species. See Appendix X for the details of model specification and evaluation.

### Model 4: Aggregate diet ordination

To visualize the net diet similarity between bumble bee species, we pooled visitation data across all samples within bumble bee species and plant species, then performed a principal coordinates analysis based on Horn-Morisita distance metric. To test the significance of tongue length class as a predictor of diet similarity, we performed a permutation test.

### Model 5: The variation of morphological bias through elevation and time

To explore the plasticity of floral biases, we first subsetted our data to include only the dominant bumble bee species (B. hortorum, B. pascuorum, B. pratorum, B. soroeensis, B. terrestris/lucorum, and B. wurflenii) and the most visited floral morphotypes (groups 3, 5, 6, and 7). We then constructed a generalized additive model (GAM) in which binary visitation was modeled as a tensor-product smooth of elevation and day-of-year, grouped by floral morphotype. This response was conditioned on linear interaction terms for bumble bee abundance : percent flower cover and percent flower cover : floral morphotype and a random intercept term for bumble bee species : floral morphotype. See Appendix X for the details of model specification and evaluation.

## Software

All analyses were conducted in `R` [@R_Core_Team_2020]. Bayesian multilevel models were implemented in Stan (Stan Development Team 2021) via the package `brms` [@Burkner_2017; @Burkner_2018], and visualizations were constructed with the packages `tidybayes` [@Kay_2020] and `ggdist` [@Kay_2021]. Ordination analysis and permutation testing were performed with the package `vegan` [@Oksanen_2019] and visualized using the packages `ggplot2` `vegan` [@Wickham_2016] and `ggvegan` (Simpson 2019). GAMs were implemented with the package `mgcv` [@Wood_2017] and visualized with the package `mgcViz` [@Fasiolo2018aa]. Data handling and visualization was performed with the `tidyverse` suite of R packages [@Fitzpatrick_2020]. Annotated R code is available in the Supplementary Material.


# RESULTS

## Summary of visitation and floral survey data

<!-- Modify this section to describe how many observations/samples were dropped due to the non-alignment of visitation and survey data -->

We recorded a total of 12,918 bumble- bee-flower interactions (excluding males) over the three years of our study. The metaweb across all sites and dates consisted of 16 bumble bee species (with *Bombus terrestric/lucorum* and *Bombus psithyrus* species groups lumped into species groups), 163 plant species (110 genera, 37 families), and 736 unique bumble- bee-plant interaction pairs. Five species --- *B. pascuorum*, *B. pratorum*, *B. soroeensis*, *B. terrestris-lucorum*, and *B. wurflenii* --- accounted for the bulk of overall bumble bee abundance and were present over the entire elevation range of our sites **(Figure 1)**. *B. hortorum*, *B. jonellus*, and *B. psithyrus* (species-group) were widespread but at lower abundance. *B. monticola* occurred at moderate abundance but was rare or absent at sites below 1000 m. The remaining species were relatively rare in our study system. *B. mendax*, *B. mucidus*, and *B. pyrenaeus* occurred primarily above 1500 m, consistent with their known affinity for high elevation habitats [@Rasmont2010-ci]. *B. gerstaeckeri*, *B. hypnorum*, *B. lapidarius* were sparsely distributed across mid-elevation sites. *B. humilis* was recorded only three times and only in one year, so we omitted it from all analyses.

Floral surveying yielded a total of 352 plant species, representing 191 genera and 52 families. Of these, 155 species from 103 genera and 35 families were observed to be visited by bumble bees. Eight species --- *Rubus idaeus*, *Rosa canina*, *Juniperus communis*, *Larix decidua*, *Salix* sp., *Caltha palustris*, *Pulmonaria officinalis*, and *Rheum barbarum* --- were recorded in visitation data but not recorded during floral surveying. Each accounted for no more than 4 visits in total over the three years of our study, and they were omitted from the analysis of floral abundance.

![(\#fig:fig1)Each floral morphotype is exemplified by a representative species common in our study system. Illustrations are modifications of public domain works obtained from www.plantillustrations.org. For details of original sources, see Appendix X.](../output/ktype_abund_gam_mod.png)

## Overall visitation probability by floral morphotype

Lip (group 5) and flag (group 6) flowers had the highest overall visitation probability, with a conditional mean visitation probabilities of 0.36 (95% CI = 0.30–0.42) and 0.33 (95% CI = 0.27–0.39), respectively (Figure 1). In second place, head flowers (group 7) had a mean conditional visitation probability of 0.19 (95% CI = 0.16–0.24). Bell (group 3) flowers were a distant third with a conditional mean visitation probability of 0.092 (95% CI = 0.07–0.12), followed by stalk-disc (group 4), funnel (group 2), and disc flowers (group 1) with conditional visitation probabilities of 0.07 (95% CI = 0.05–0.09), 0.05 (95% CI = 0.04–0.07), and 0.05 (95% CI = 0.03–0.07), respectively. Unclassified (group 0) flowers, consisting mainly of wind-pollinated grasses, were the least-visited, with a conditional mean visitation probability of 0.03 (95% CI = 0.02–0.04). Our model explained approximately 35% of the total variance (97.5% CI = 34–36%).


![(\#fig:fig2)Posterior conditional mean visitation probability by floral morphotype. Point-interval plots depict 66\% and 95\% CIs.](../output/brm_pvis01_plot1.pdf)


## Visitation probability by flower morphotype, bumble bee species, and bumble bee tongue length class

Descriptions of tongue classes are based on Model 2 and refer to Figure 2, and descriptions of species are based on Model 3 and refer to Figure 3. Both models explained approximately 28% of total variance (97.5% CI = 27–29%). As in Model 1, visitation rates refer to conditional visitation probabilities.

Disc (group 1) flowers were visited at similar rates across tongue length classes. At the species level, disc flowers were visited at relatively high rates by the short-tongued B. terrestris/lucorum and the medium-tongued Bombus monticola. Visitation to funnel (group 2) flowers was dominated by long-tongued species, an effect driven entirely by the Aconitum (group 2) specialist B. gerstaeckeri. Bell (group 3) flowers were visited at a higher rate by short-tongued species than by long-tongued species, with medium-tongued species exhibiting intermediate visitation rates. At the species-level, bell flowers were visited at highest probability by the short-tongued species B. jonellus, B. pyrenaeus, and B. soroeensis, and by the medium-tongued species B. hypnorum and B. monticola. Stalk-disc (group 4) flowers were visited at similar rates across all tongue length classes, and there was also broad overlap in visitation rates at the species level. Lip (group 5) flowers were also visited at similar rates across tongue length classes, but at the species level B. pascuorum exhibited an uniquely high preference for lip flowers. Flag (group 6) flowers were visited at a higher rate by medium-tongued species than by short-tongued species, with long-tongued species exhibiting intermediate visitation rates. Flag flowers generated strong separation in visitation rate at the species level, receiving very high visitation from the medium-tongued B. pascuorum and B. mucidus, moderately high visitation from the short-tongued B. terrestris/lucorum and B. lapidarius, the medium-tongued B. monticola, and the long-tongued B. mendax and B. hortorum, and low visitation by the short-tongued B. soroeensis, B. pyrenaeus, B. psithyrus, B. pratorum, and B. jonellus, the medium-tongued B. hypnorum, and the long-tongued B. gerstaeckeri.


![(\#fig:fig3)Posterior conditional mean visitation probability by Kugler morphotype and bumble bee tongue length class. Point-interval plots depict 66\% and 95\% CIs. Letters indicate group differences based on evidence ratio at alpha = 0.05. Note that the range of the x-axis has been rescaled for each panel to emphasize the comparison of visitation probabilities across tongue length classes.](../output/brm_pvis00_plot1.pdf)

![(\#fig:fig4)Posterior conditional mean visitation probability by Kugler morphotype and bumble bee species. The tongue length class of each species is color-coded. Point-interval plots depict 66\% and 95\% CIs. Note that the range of the x-axis has been rescaled for each panel to emphasize the comparison of visitation probabilities across bumble bee species](../output/brm_pvis02_plot1.pdf)

## Diet ordination

Axis 1 accounted for 24% of total diet variance and primarily separated long- and medium-tongued species from short-tongued species (Figure 4). Long- and medium-tongued species scored low on Axis one and were associated with the lip flowers (group 5) Ajuga reptans, Rhinanthus glacialis, and Stachys alopecuros) and flag flowers (group 6) Anthyllis vulneraria, Lotus corniculatus, and Trifolium pratense. Short tongued species scored high on Axis one and were associated with the bell flowers (group 3) Campanula scheuchzeri, Erica carnea, and Geum rivale and the head flowers (group 7) Carduus defloratus, Centaurea jacea, and Phyteuma orbiculare. Axis 2 accounted for 19% of total variance and primarily explained variation within short-tongued species. Short-tongued species scoring high on Axis 2 were associated with Campanula scheuchzeri while those scoring low on Axis 2 were associated with Erica carnea. Short-tongued species scoring close to zero on Axis 2 also scored lower on Axis 1 (closer to medium- and long-tongued species) and were associated with head flowers. Tongue length class explained 22% of diet variance (Adonis permutation test, F = 1.68, p = 0.006). While short-tongued bees appeared most broadly spread in ordination space, there were no significant differences in dispersion across tongue length classes (F = 0.62, p = 0.55).

![(\#fig:fig5) Principal coordinates ordination of bumble bee visitation pooled across all samples.](../output/pcoa_aggregate.png)


## The variation of morphological bias through elevation and time

Model 5 explained 33.8% of total variance (adjusted R2 = 0.37). All tensor product smooths (Figure 5) were significant at p < 0.001 except for the one representing B. pratorum and flag (group 6), which was more weakly significant at p = 0.02. As in Models 1-3, the results of Model 5 should be interpreted in terms of conditional visitation probability.

Bell (group 3) flowers were visited at highest rates early in the season and at low elevation by all species except B. hortorum, whose visitation rate to bell flowers peaked in mid-season and at mid-elevation. Notably, the medium-tongued B. pascuorum and B. wurflenii also exhibited relatively high visitation to bell flowers at mid-season and mid-elevation, while the visitation rates of the short-tongued B. pratorum, B. soroeensis. and B. terrestris lucorum were more concentrated at the early-season and low-elevation peak.

Lip (group 5) flowers exhibited roughly the mirror image of the pattern of visitation rate seen in bell flowers, with peak visitation late in the season and at low elevation, but relatively high visitation extending toward mid-season and mid-elevation. This pattern was consistent for all bumble bee species.

The pattern of visitation rates to flag (group 6) flowers fell generally along the lines of tongue length class. The short-tongued B. pratorum and B. soroeensis exhibited peak visitation early in the season and at low elevation. The long-tongued B. hortorum and the medium-tongued B. pascuorum and B. wurflenii exhibited peak visitation at mid-season, though B. pascuorum exhibited a bimodal pattern in elevation, with peaks around 1000 m and 1750 m, while the visitation rate of B. hortorum declined monotonically from high elevation to low. The short-tongued B. terrestis/lucorum broke the tongue-length pattern, exhibiting a diffuse pattern of visitation rate more similar to that of the long- and medium-tongued species than the other short-tongued species.

Visitation rates for head (group 7) flowers were more or less bimodal for all species, with a first peak in mid-season and low elevation and a second peak in late-season and high elevation.

![(\#fig:fig6)GAM analysis of bias as a function of the tensor-product interaction of elevation and day of year, parsed by bumble bee species and floral morphotype and conditioned on the abundance of each bumble bee species and the relative abundance of each floral morphotype. The color ramp of each smooth has been scaled to equal color range, thereby obscuring the absolute magnitude of effect. To avoid cluttering the figure, we have opted not to include a unique color ramp key for each plot, since the focus of our inference is on relative rather than absolute patterns.](../output/ktype_lability_gam_mod.png)

# DISCUSSION

## Foraging bias and floral morphology

Our study is not the first to compare bumble bee visitation rates across Kugler’s (1970) floral morphotypes (see Neumayer and Paulus 1999; Schneller et al. 2014), but to our knowledge it is the first to enable the inference of bias across morphotypes by controlling statistically for the neutral effects of abundance.

When bumble bee visitation is pooled across all species, it is clear that lip, flag, and head flowers dominate the floral community in terms of conditional visitation probability. Patterns of pooled visitation are, of course, driven disproportionately by the behavior of the most abundant bumble bee species, which in our system were B. hortorum, B. pascuorum, B. pratorum, B. psithyrus, B. soroeensis, B. terrestris/lucorum, and B. wurflenii. Visitation to lip and flag flowers was driven especially by the medium-tongued B. pascuorum, the most abundant bumble bee in our study system, and secondarily by the medium-tongued B. wurflenii. Visitation to head flowers was driven mainly by the short-tongued B. soroeensis and B. psithyrus.

The influence of tongue length on visitation bias was pronounced for some morphotypes but not for others. In the cases of disc, and stalk-disc, and unclassified flowers, visitation rates were uniformly low across all tongue length classes. The opposite was true of lip flowers, which were visited at similarly high rates by all tongue length classes. Bell flowers and head flowers were visited chiefly by short-tongued bees, consistent with the shallow corollas of these morphotypes. In contrast, funnel flowers, with their deep corollas were visited chiefly by long-tongued bees, though this effect was driven disproportionately by B. gerstaeckeri, which specializes on the funnel flowers of the genus Aconitum. Flag flowers were visited at the highest rate by medium-tongued species, and especially by B. pascuorum and B. mucidus.

While the variation in floral bias across tongue length classes we observed largely comports with the established relationship between tongue length and corolla depth, variation in bias within tongue length classes demonstrates the insufficiency of simple morphological matching to explain bumble bee foraging patterns. Nectar robbing is one mechanism by which the relationship between tongue length and corolla depth might be relaxed, but its influence in our data is unclear. Of the bumble bees in our study system, B. terrestris/lucorum and especially B. wurflenii are known to engage frequently in nectar robbing (Brian 1954). In the case of B. wurflenii, however, its morphotype affinities did not differ markedly from those of other medium-tongued species, and it may be that the condition of medium-tongue length alone confers enough to foraging flexibility in our study system to dampen the divergent effects of nectar robbing. Nectar robbing might, however, account for the deviation of B. terrestris/lucorum from most other short-tongued species with respect to flag flowers. With the exception of B. lapidarius, all short-tongued species besides B. terrestris/lucorum eschewed flag flowers, while B. terrestris/lucorum visited them at rates comparable to medium- and long-tongued species. Indeed, B. terrestris/lucorum was, on the whole, an idiosyncratic forager, demonstrating a relatively high affinity for disc flowers and unclassified flowers (particularly Plantago lanceolata) that were largely ignored by other species. Perhaps the most striking divergence of floral bias within tongue length class, though, was the divisive effect of bell flowers on short-tongued bumble bees. B. pyrenaeus, B. soroeensis, B. pratorum, and B. jonellus visited bell flowers at high rates, while B. terrestris/lucorum, B. psithyrus, and B. lapidarius visited them relatively rarely. A similar split among the short-tongued species occurred with respect to head flowers, which were visited at high rates by B. soroeensis and B. psithyrus but at relatively low rates by other short-tongued species. There are no obvious morphological distinctions within short-tongued bumble bees that would account for these divergent affinities; ecological drivers such as competition or differences in spatial or temporal activity patterns are more likely explanations, as will be revisited below in our discussion of overall diet similarity.


## Diet similarity

Our finding that bumble bee diet ordination produces a primary axis corresponding to tongue length — and associated with affinity for legumes — and a secondary axis describing the dispersion of short-tongued species is consistent with studies in Poland (Dave Goulson, Lye, and Darvill 2008), France (Durieux 2000), and England (D. Goulson et al. 2005) that analyzed bumble bee diet by similar methods, suggesting that this may be a general pattern among bumble bee communities, at least within temperate Europe. While our test of beta-dispersion did not reveal significant differences between tongue length classes, the spread of short-tongued bees along the secondary axis of our ordination, together with their divergence with respect to visitation rates on bell and head flowers described above, suggests that resource partitioning may be more pronounced among short-tongued bumble bees than among their medium- and long-tongued counterparts. It is tempting to invoke interspecific competition as a cause of the apparent niche partitioning among short-tongued bumble bees, in agreement with previous studies that have demonstrated the potential for severe competition between bumble bees of similar tongue length (e.g. Heinrich 1976; Inouye 1978; Pyke 1982). But the importance of competition in bumble bee communities is not universal (e.g. Ranta and Vepsäläinen 1981; Williams 1989), and in the rich wildflower meadows of our (sub)alpine system, bumble bee populations are more likely limited by harsh abiotic conditions than by food availablility. Alternatively, it may be that short-tongued bees in our study system differed from one another in their affinities for specific habitats (Schneller et al. 2014) or spatial configurations of resource patches (Sowig 1989), and their divergent morphological biases arose secondarily from these effects.


## The lability of bias through elevation and time

While the intraspecific lability of bumble bee foraging biases has been demonstrated under laboratory conditions (Ings, Raine, and Chittka 2009; Maharaj et al. 2019) and in a localized field study (Raine and Chittka 2007), ours is the first study, to our knowledge, to characterize the in situ lability of foraging biases within multiple species through elevational and temporal gradients of floristic turnover. In our analysis, it is evident that morphological biases in bumble bees, insofar as they are captured by the morphotypology used in our study, are not fixed but rather tuned to floristic (and possibly abiotic) context and therefore dynamic through elevation and time. Notably, bumble bee species differed from one another not just in their overall bias toward given morphotypes but also in when and where they were biased toward a given morphotype. This is especially evident in the case of flag flowers, for which positive bias was maximized in mid-season for B. hortorum, B. pascuorum, B. terrestris/lucorum, and B. wurflenii but in early-season for B. pratorum and B. soroeensis. This pattern corresponds to the overall biases of these species toward flag flowers, which were primary resources for B. hortorum, B. pascuorum, B. terrestris/lucorum, and B. wurflenii but only marginal resources for B. pratorum and B. soroeensis.


# CONCLUSIONS

For bumble bees, floral bias is a labile product of innate preferences, the soft constraints of morphological compatibility, and the dynamic ecological context in which foraging occurs. By conditioning bumble bee visitation data upon concurrent floral abundance data across spatially and temporally replicated samples, we demonstrate the strength bumble bee floral bias, its functional underpinnings, its variation across species, and its lability within species in response to environmental variation. Such tuning of floral bias has obvious adaptive significance for bumble bees, which are characteristically — and most likely primordially — associated with mountain habitats (Williams 1985; Hines 2008), where steep temperature gradients and short growing seasons generate extreme spatial and temporal floral turnover.
   



# ACKNOWLEDGEMENTS

We thank D.J. McNeil for helpful conversations during data analysis.


# REFERENCES

<!-- ```{r write_citations, cache=FALSE, include=FALSE} -->
<!-- write.bibtex(file = "knitcitations.bib") -->
<!-- ``` -->

<div id ="refs"></div>
