---
title: "Longitudinal Analysis of Child Microbiomes and Brain Development"
author: "Kevin Bonham, PhD"
institute: "Wellesley College"
topic: "Gaussian Process and Longitudinal Microbiome Analysis"
theme: "Luebeck"
fonttheme: "professionalfonts"
mainfont: "Hack Nerd Font"
fontsize: 10pt
urlcolor: red
linkstyle: bold
# titlegraphic: notes/assets/Alexa_brain.png
# logo: img/aleph0-small.png
date: "2021-12-17"
lang: en-US
section-titles: false
toc: true
---


# Background

## The microbiome develops rapidly in the first year of life

::: columns

:::: column
- No (or little) prenatal microbiome
- Perinatal microbiome strongly influenced by delivery,
  antibiotics, diet (breast milk vs formula)
- Gut microbiome is seeded by delivery (birth canal),
  skin of caregivers, and milk source
- Slowly increasing diversity, conversion to "adult-like" microbiome
  after transition to solid food
::::

:::: column

![](https://i.imgur.com/PgA5D73.png)

::::

:::

## The brain develops rapidly in the first years of life

::: columns

:::: column

- By age 5, brain has reached 85% of adult size
- Gross pattens of axonal connections are established,
  and near-adult levels of myelenation also achieved by 5 years of age
- Development driven in part by environmental exposures
  (care-giver attention, stress, diet, etc)

::::

:::: column
![](notes/assets/Brain-development.jpg)
leelanauearlychildhood.org/brain-development
::::

:::


## The gut-brain-microbiome axis: bi-directional interactions between microbiome and nervous system

![](notes/assets/gut-brain-axis.jpg){ height=80% }
https://doi.org/10.1016/j.tins.2013.01.005

## Evidence of microbiome-brain interactions from autism in humans

![](https://i.imgur.com/IQ2ZLEG.png){ height=75% }
https://doi.org/10.1038/s41598-018-36430-z


## Evidence of microbiome-brain interactions from autism in humans

![](https://i.imgur.com/mtmGxrq.png){ height=75% }
https://doi.org/10.1038/s41598-018-36430-z

## Evidence of microbiome-brain interactions from "autism" in mice

![](https://i.imgur.com/ztHCf8L.png){ height=90% }
http://dx.doi.org/10.1016/j.cell.2013.11.024

## Open questions

###

What effects does the microbiome on have on normal brain development?

. . .

###

What are the links between microbial metabolism,
  gut metabolites, and neurocognitive function?

# The RESONANCE cohort of child brain development

## Cohort design

::: columns

:::: {.column width=40%}

- Shotgun metagenomic sequencing
- LCMS Metabolomics
- Neuroimaging (MRI)
- Cognitive assessments (IQ-like)
- Sleep, genetics, nutrition
- Other clinical covariates

:::

:::: {.column width=60%}

![](notes/assets/microbiome_cohort.png){ height=90% }

::::

:::


##

### Starting with cross-sectional data:

Are specific taxa, gene functions, or metabolites associated with neurocognition?

## Some taxa are associated with cognition

![](notes/assets/cogscores.png){ height=80% }
Upper and lower quartiles of cognitive function

## Potentially neuroactive genes are associated with cognition and social responsiveness

::: columns

:::: column

![](notes/assets/cogscore-fsea.png){ height=90% }
FSEA of neuroactive genes for cognitive function

::::

. . .

:::: column

![](https://i.imgur.com/xlkoV18.png){ height=90%}
FSEA of neuroactive genes with SRS2 scores

in collaboration with Hannah Laue
::::

:::

Valles-Colomer _et. al._ (2019)

## Potentially neuroactive genes are associated with cognition and social responsiveness

::: columns

:::: column

![](notes/assets/cogscore-fsea-box.png){ height=90% }
FSEA of neuroactive genes for cognitive function

::::

:::: column

![](https://i.imgur.com/xlkoV18.png){ height=90%}
FSEA of neuroactive genes with SRS2 scores

in collaboration with Hannah Laue
::::

:::

Valles-Colomer _et. al._ (2019)

## Infants have very different metabolomes than their mothers

::: columns

:::: column

![](notes/assets/metabolites_pcoa.png){ height=80% }
PCoA of metabolites measured with 4 methods

::::

. . .

:::: column

![](notes/assets/gaba-glutamate.png){ height=80% }

::::

:::

##

### Question:

How are microbial metabolic potential and gut metabolome linked,
esp for neuroactive genes?

. . .

### GABA and Glutamate are particularly important in neuronal development in infants

And there are known links between these molecules in the gut and brain function

## Microbial metabolism of GABA and Glutamate metabolism is not significantly associated with gut concentrations

::: columns

:::: column

![](notes/assets/gaba_genes_metabolites.png){ height=80% }

::::

. . .

:::: column

![](notes/assets/glutamate_genes_metabolites.png){ height=80% }

::::

:::

## Relationships between microbial metabolism and molecules

![](https://i.imgur.com/uhCgQBI.png)

## Limitations due to pandemic

![](https://i.imgur.com/7aJk3zC.png)

## Conclusions

- The early childhood microbiome is associated with brain development
- Possible effects through SCFAs and metabolism of other neuroactive molecules
- Don't try to do longitudinal human cohort stuff during a global pandemic

### 

"When self and other are the same, mind and ~~dharmas~~ microbes are
one."

~ Hongzhi Zhengjue _Cultivating the Empty Field_



# The GP tool - GaPLAC

## Which would you choose?

![](notes/assets/model_compare-anim1.png){ height=90% }

## Which would you choose?

![](notes/assets/model_compare-anim2.png){ height=90% }


<!-- ## Linear models identify parameters that satisfy linear relationship of random variables

![](notes/assets/LMM_anim-1.png){ height=90% }
Fitting a linear model

## Linear models identify parameters that satisfy linear relationship of random variables

![](notes/assets/LMM_anim-2.png){ height=90% }
Fitting a linear model

## Linear models identify parameters that satisfy linear relationship of random variables

![](notes/assets/LMM_anim-3.png){ height=90% }
Fitting a linear model

## Linear models identify parameters that satisfy linear relationship of random variables

![](notes/assets/LMM_anim-4.png){ height=90% }
Fitting a linear model -->

<!-- ## GPs identify functions that satisfy covariance relationship between points

![](notes/assets/GP_anim-1.png){ height=90% }
Fitting a Gaussian process model

## GPs identify functions that satisfy covariance relationship between points

![](notes/assets/GP_anim-2.png){ height=90% }
Fitting a Gaussian process model

## GPs identify functions that satisfy covariance relationship between points

![](notes/assets/GP_anim-3.png){ height=90% }
Fitting a Gaussian process model

## GPs identify functions that satisfy covariance relationship between points

![](notes/assets/GP_anim-4.png){ height=90% }
Fitting a Gaussian process model -->

## Gaussian Process models (GPs)

- Collection of arbitrarily many random variables linked by a covariance function (kernel).

. . .

- Entirely described by this kernel, which dictates how “nearby”
  points relate to each other.

. . .

- Typically over time, but can be any other continuous variable).

. . .

- More constrained that freeform nonparametric Bayes, more flexible than regression.

. . .

![](https://i.imgur.com/wz9wOel.png){ height=40% }
source: Wikipedia

## Linear models vs GPs

::: columns

:::: column

![](notes/assets/LMM_anim-4.png){ height=90% }
Fitting a linear model 

::::

:::: column

::::

:::

## Linear models vs GPs

::: columns

:::: column

![](notes/assets/LMM_anim-4.png)
Fitting a linear model 

::::

:::: column

![](notes/assets/GP_anim-4.png)
Fitting a GP
::::

:::

## GPs in HMP1-II to model different sources of variation

![](https://i.imgur.com/rkOoezk.png){ height=70% }
Credit: Jason Lloyd-Price https://doi.org/10.1038/nature24485

## GPs in HMP1-II to model different sources of variation

![](https://i.imgur.com/tWNlOzP.png){ height=70% }
Credit: Jason Lloyd-Price https://doi.org/10.1038/nature24485

## GPs in the literature to model other dynamic processes

![](notes/assets/gps-literature.png){ height=80% }

## GPs could provide more realistic model for longitudinally sampled microbiomes

![](notes/assets/model_compare-anim2.png){ height=90% }

## Could GPs better handle transitions between infant and adolescent communities?

![](notes/assets/age_taxon_ratios.png){ height=75% }


## (Ga)ussian (P)rocess models for (L)ongitudinal (A)nalysis of (C)ommunities

![](https://i.imgur.com/bEapMwf.png){ height=50% }

https://en.wikipedia.org/wiki/Geplak

. . .

### Goal

Make fitting GPs as easy as fitting LMs.

## Current Functionality - Sample (simulate)

::: columns

:::: column

![](https://i.imgur.com/9hIyKNV.png){ width=90% }

![](notes/assets/sqexpplot.png){ width=90% }

::::

. . .

:::: column

![](https://i.imgur.com/RLxwIAy.png){ width=90% }

![](notes/assets/ouplot.png){ width=90% }

::::

:::

## Current Functionality - Sample (simulate)

::: columns

:::: column

![](https://i.imgur.com/oqVikKB.png){ width=90% }

![](notes/assets/linearplot.png){ width90=% }

::::

:::: column

::::

:::

## Current Functionality - MCMC (fit / infer parameters)

![](https://i.imgur.com/mhxp9uj.png)

. . .

![](https://i.imgur.com/s7yt5RT.png)

## Current Functionality - Select (compare models)

![](https://i.imgur.com/V6NAfUb.png)

. . .

![](https://i.imgur.com/mvmnMOk.png){ height=40% }

## Performance on simulated data - it works!

::: columns

:::: column

**Synthetic Data:**

- 200 subjects, 300 “bugs”.
- One bug “_Faecalibacterium prausnitzii_”
  linearly associated with continuous variable "diet."
- Modeled from real HMP gut microbiomes using SparseDOSSA 2.
- Diet = simple continuous variable.
- Three (randomly sampled) time points / individual.
- Compared Bayes factor (essentially, fit) for model
  with or without diet variable
::::

:::: column

![](https://i.imgur.com/WgV3Mnw.png){ height=80%}

Special thanks: Meg

::::

:::

## Conclusions from real-world data


![](https://i.imgur.com/jCi4kpW.png)

Special thanks: Tobyn

## Conclusions from real-world data

![](https://i.imgur.com/GgGZV8B.png)

- Comparing model: `Bug ~ SqExp(time) * Categorical(subject)`
- With or without `Linear(diet)`

Special thanks: Tobyn

## Next steps

- Figure out what's going on with MCMC sampling
- Incorporate non-Gaussian likelihoods (similar to link functions in LMs)
- Allow different AD backends / sampling methodology (eg NUTS)
- Benchmark against MaAsLin on simulated data
- Compare against MaAsLin on real-world data

## Thanks!

::: columns

:::: column

### Vanja Klepac-Ceraj

Shelley McCann

Sophie Rowland

_Danielle Peterson_

_Lauren Tso_

Anika Luo

Annelle Abatoni

Alexa Gross

::::

:::: column

### Juliette Madan

_Hannah Laue_

### Curtis Huttenhower

_Jason Lloyd-Price_

_Tobyn Branck_

_Meg Short_

Andrew Ghazi

Eric Franzosa

::::

:::

http://www.infinitecuriosity.org/vizgp/