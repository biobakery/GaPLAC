
# Guide

This guide is intended to provide an overview of the basic workflow using GaPLAC. A complete command reference is provided below.

## Installation

## Generating some sample data

***generate and plot some sample data using `sample`***

## Fitting parameters

***run mcmc to recover parameters***

## Comparing models

***run select to determine which model we should use***

## Making predictions

## Diagnostic plots

# Command reference

## `sample`

## `mcmc`

## `select`

## `predict`

## `fitplot`

# Input and output

## Input variables

The variables which can be used in the formulas are provided to GaPLAC using the `--data` parameter, in the form of a set of tab- or comma-separated tables. The column/row names of these tables are used as the variable names.

***non-conforming variable names***.

***--data format***

## Target input variables

***--atdata and --at formats***

## Output

***--output format***

***--plot format***

# Gaussian process formulas

Formulas specify the form of the covariance function of the Gaussian Process, its data likelihood, as well as the values of any hyperparameters and their priors. Because of this, the formula effectively encodes the hypothesis you are testing when using GaPLAC.

Formulas have the following syntax:

```
observation [: likelihood] ~| covariance_function
```

Formulas are designed to resemble model formulas in R to allow those familiar with R's model formulas to gain some intuition into working with these formulas. However, it is important to remember that model formulas in R are used to specify fixed effects, while formulas here model covariance. To distinguish these, note that formulas here include `~|` rather than simply `~`, following the usual specification of random effects in R's formulas using `|`.

The `observation` can be any Julia expression of the input variables. The `covariance_function` is an expression of covariance function components, and has more restricted syntax. It is recommended that you test what types of functions are described by a covariance function before using it in your workflow using the `sample` command.

## Basic covariance functions

### Constant

### Categorical

### Ornstein-Uhlenbeck

### Squared-exponential

`gaplac sample "y ~| Cat(person) * SExp(time; l=1)" --at "person=1:4;time/person=-5:0.1:5" --plot sample_plot.png --plotx time:person`

### Periodic

### Linear




## More complex examples

`gaplac sample "y ~| Cat(person) * <covariance function>" --at person=1:4;time=-5:0.1:5 --plot sample_plot.png --plotx time:person`


## Variance parameters

***discuss how variance parameters are handled***

## Data likelihoods

### Predictions

![picture](img/abc.png)

