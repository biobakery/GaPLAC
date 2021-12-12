---
title: "Longitudinal Analysis of Child Microbiomes and Brain Development"
author: "Kevin Bonham, PhD"
institute: "Wellesley College"
topic: "Gaussian Process and Longitudinal Microbiome Analysis"
theme: "Frankfurt"
colortheme: "beaver"
fonttheme: "professionalfonts"
mainfont: "Hack Nerd Font"
fontsize: 12pt
urlcolor: red
linkstyle: bold
aspectratio: 169
# titlegraphic: img/aleph0.png
# logo: img/aleph0-small.png
date:
lang: en-US
section-titles: false
toc: true
---

# Background

## The microbiome develops rapidly in the first year of life

- No (or little) prenatal microbiome
- Perinatal microbiome strongly influenced by delivery,
  antibiotics, diet (breast milk vs formula)
- Gut microbiome is seeded by delivery (birth canal),
  skin of caregivers, and milk source
- Slowly increasing diversity, conversion to "adult-like" microbiome
  after transition to solid food

## The brain develops rapidly in the first years of life

- By age 5, brain has reached 85% of adult size
- Gross pattens of axonal connections are established,
  and near-adult levels of myelenation also achieved by 5 years of age
- Development driven in part by environmental exposures
  (care-giver attention, stress, diet, etc)

## The gut-brain-microbiome axis

<!-- Some images from papers -->

## Evidence of microbiome from pathology

- esp autism
- Autism mouse models
- other behavioral outputs

## Open questions

- effects of microbiome on normal brain development
- links between microbial metabolism, gut metabolites, and neurocognitive function

# The RESONANCE cohort of child brain development

## Cohort design

<!-- model image -->

## Data Collection

- shotgun metagenomic seuqencing
- neuroimaging
- cognitive assessments (IQ-like)
- sleep
- genetics
- Other clinical covariates

## Early evidence for effects

<!-- cog score data -->

## Functional analysis (FSEA)

<!-- cog score data -->

## Limitations due to pandemic

<!-- Sample collection -->

# Methods for longitudinal analysis

## Background on linear models

## Application of linear models to longitudinal data (mixed effects models)

## Challenges with microbiome data

- unevenly sampled timepoints
- non-uniform sparsity (eg kids)

# The GP tool - GaPLAC

## Background on GPs

- HMPI-II
- other gp papers

## Motivation for using GPs

- Example: irregular longitudinal sampling
- Example: Bugs that disappear or appear with age

## (Ga)ussian (P)rocess for (L)ongitudinal (A)nalysis of (C)ommunities

- **Goal**: Make fitting GPs as easy as fitting LMs

<!-- Example of MixedEffects Model math vs R code -->

## Model definition



## Current Functionality - Sample

<!-- Example of sample -->

## Current Functionality - MCMC

<!-- Example of mcmc -->

## Current Functionality - Select

<!-- Example of select -->

## Performance on simulated data



## An H2

Some other stuff

## Do columns work out of the box?

::: columns

:::: column

Left column text.

Another text line.

::::

:::: column

- Item 1.
- Item 2.
- Item 3.

::::

:::