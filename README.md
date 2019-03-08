# Guide

This guide is intended to provide an overview of the basic workflow using GaPLAC. A complete command reference, as well as available covariance and likelihood functions are provided below.

## Installation

1. Install [Julia](https://julialang.org/).

2. Download GaPLAC's repository and unpack it somewhere.

3. Open a console in GaPLAC's root folder and run `julia`.

4. At the Julia prompt, enter `]activate .`, then `instantiate`. This will install the required packages. Use backspace to get back to the normal Julia prompt, and run `exit()` to quit Julia.

5. If on Mac/UNIX, to use `./gaplac ...` format, you may need to run `chmod u+x ./gaplac`.

## Generating some sample data

GaPLAC has five main commands to work with: `sample`, `mcmc`, `select`, `predict`, and `fitplot`. We will first look at `sample`, which simply draws a sample from a Gaussian Process. This can be helpful to visualize the kinds of functions described by a particular GP, or in this case to provide some sample data to work with in the later functions.

Run the following command from the GaPLAC root folder:

```
./gaplac sample "y :~| SExp(x; l=1)" --at "x=-5:0.1:5" --plot gp_sample.png
```

This will produce a large amount of output to the console, and should also produce a plot in `gp_sample.png` which resembles:

![Wavy line](img/guide1.png)

Let's look at each of the pieces of the command:

- `"y :~| SExp(x; l=1)"`: This is the GP formula, much like a model formula in R. In this case, the output (`y`) is modeled as a GP with a Squared-Exponential covariance function (`SExp`) with a lengthscale (`l`) of `1`. Note also the `:` in `:~|`. Normally, a data likelihood (see below) can be specified between the `:` and the `~`, but here we don't specify anything, and the GP will be modeled _without_ a likelihood. This effectively allows us here to draw samples from _just_ the Gaussian Process described by the formula, without additional variation from the likelihood.
- `--at "x=-5:0.1:5"`: This tells GaPLAC what values of `x` to sample the GP at.
- `--plot gp_sample.png`: Output the pretty plot here.

Try changing the lengthscale of the `SExp` term. How does this affect the function? Try adding other components (the full list is at the end of this document) by adding them to the formula, such as an Ornstein-Uhlenbeck process (`OU(x; l=1)`), or simply some `Noise`.

Now let's generate a smaller set of data at some randomly chosen `x` coordinates, and store the results in a file instead of printing to stdout:

```
./gaplac sample "y :~| SExp(x; l=1.5)" --at "i=1:50;x=Uniform(-5,5)" --plot gp_sample.png --output data.tsv
```

Look at the contents of `data.tsv`. It should contain two columns: `x` and `y`, and the rows are not sorted in any way. We will use this data for the next command.


## Fitting parameters

We are usually interested in the parameters of the covariance function which best fits some data. This is accomplished with the `mcmc` command in GaPLAC, which will produce a [MCMC chain](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) of samples from the posterior distribution of the model parameters.

Try running the following command:

```
./gaplac mcmc "y ~| SExp(x)" --data data.tsv --output mcmc.tsv --samples 500
```

First, note that the model formula omits the additional `:` before the `~`. This will therefore default to a Gaussian likelihood. During inference, you will almost always want to have at least some kind of noise/error term which allows for some additional variation from the GP, as fitting without this term can lead to computational difficulties.

Now examine the output file `mcmc.tsv`. Instead of showing a relationship between `x` and `y`, as the previous run of `sample` did, this file should contain several columns of parameter values, as well as the all-important final column containing the log of the unnormalized posterior density for the sample (think of this something like the goodness of fit). In particular, take a look at the covariance function parameter `l`, which we set above in the model formula to `1.5`, but which we did not tell the `mcmc` command. If all worked well, the mean of this parameter should converge to, and hover around the true value of `1.5`.

## Making predictions

Armed with some likely parameter values, we can now use the trained model to make some predictions for unobserved values. Try running the following:

```
./gaplac predict "y ~| SExp(x)" --data data.tsv --mcmc mcmc.tsv --at "x=-8:0.1:8" --output prediction.tsv
```

This will produce a prediction for the values of the process for `x` values from -8 to 8. Since the dataset ranges from -5 to 5, some of these are far from actual data points. Try plotting the various quantiles in `prediction.tsv` (the `Q` columns) and observe how the prediction relates to the presence of data points.

## Comparing models

`mcmc` fits a single model, but how do we know it's the right model? Maybe another model would be a better fit? To answer this question, we can run `mcmc` again with the other model. In this case, let's test a different kind of time-varying process called an Ornstein-Uhlenbeck (OU) process:

```
./gaplac mcmc "y ~| OU(x)" --data data.tsv --output mcmc_ou.tsv --samples 500
```

This will give us a second set of model fit results in a new `mcmc_ou.tsv` file. Now we can use those goodness-of-fit values in each of the files to determine which of the models we believe (hint: we generated data from a Squared-Exponential covariance function, and thus we expect the OU process will perform worse. Let's test that with the `select` command:

```
./gaplac select --mcmc mcmc.tsv --mcmc2 mcmc_ou.tsv
```

This will compare the log posterior values stored in each of the MCMC chains, and summarize them as a [Bayes Factor](https://en.wikipedia.org/wiki/Bayes_factor), which is reported in log2 scale. Here, log2 Bayes Factors greater than 1 indicate that the first model (in this case the Squared-Exponential) should be preferred, while negative numbers indicate the opposite - that the second model (passed by `--mcmc2`) should be preferred.

## Diagnostic plots

GaPLAC also has the capability to automatically generate several plots showing the correspondence between the data and the model fit. These are automatically generated based on the model formula. These plots are generated with the `fitplot` command, passing it the formula, the data, and the MCMC chain:

```
./gaplac fitplot "y ~| SExp(x)" --data data.tsv --mcmc mcmc.tsv --output fitplots.pdf
```

Take a look at the output in `fitplots.pdf`.


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

### Categorical

Syntax: `Cat(x)`

> k(i,j) = 1 if x[i] == x[j], 0 otherwise

Produces a covariance of `1` between samples with the same `x`, otherwise `0`.

`./gaplac sample "y :~| Cat(k)" --at "k=1:4;x/k=-5:0.1:5" --plot sample_plot.png --plotx x:k`

Categorical covariance functions are frequently multiplied by another function from the set below. This allows the covariance structure within each category to be described by the function being multiplied in. For example, if making measurements from individual people, we would expect that the measured value from any given person does not change very rapidly (e.g. if modeling weight over time), but different people are independent (i.e. my weight does not affect your weight).

### Constant

Syntax: `Constant(σ2=magnitude)`

> k(i,j) = 1

A constant covariance with the given magnitude. The shorthand `1` can be used for this, as in `y ~| 1`.

`./gaplac sample "y :~| Cat(k) * Constant" --at "k=1:4;x/k=-5:0.1:5" --plot sample_plot.png --plotx x:k`

### Noise

Syntax: `Noise`

> k(i,j) = 1 if i == j, otherwise 0

`./gaplac sample "y :~| Noise" --at "x=-5:0.1:5" --plot sample_plot.png --plotx x`

Uncorrelated variation - adds a covariance of `1` only from a sample to itself.

### Ornstein-Uhlenbeck

Syntax: `OU(x; l=lengthscale)`

> k(i,j) = exp(abs(x[i] - x[i])/l)

`./gaplac sample "y :~| Cat(k) * OU(x; l=1)" --at "k=1:4;x/k=-5:0.1:5" --plot sample_plot.png --plotx x:k`

### Squared-exponential

Syntax: `SExp(x; l=lengthscale)`

> k(i,j) = exp(.5*(x[i] - x[i])^2/l^2)

A squared

`./gaplac sample "y :~| Cat(k) * SExp(x; l=1)" --at "k=1:4;x/k=-5:0.1:5" --plot sample_plot.png --plotx x:k`

### Periodic

Syntax: `Periodic(x; l=lengthscale, p=period)`

> k(i,j) = exp(.5*(x[i] - x[i])^2/l^2)

`./gaplac sample "y :~| Cat(k) * Periodic(x; l=1, period=2)" --at "k=1:4;x/k=-5:0.1:5" --plot sample_plot.png --plotx x:k`

Decaying periodicity can be encoded by mixing in an Ornstein-Uhlenbeck or Squared Exponential component:

`./gaplac sample "y :~| Periodic(x; l=0.3, period=1) * SExp(x; l=2)" --at "x=-10:0.1:10" --plot sample_plot.png --plotx x`


### Linear

Syntax: `Linear(x)`

> k(i,j) = x[i] * x[j]

`./gaplac sample "y :~| Cat(k) * Linear(x)" --at "k=1:4;x/k=-5:0.1:5" --plot sample_plot.png --plotx x:k`

## More complex functions

Covariance functions can be added together, in which case the functions they represent are added together, or they can be multiplied together, in which case their effects attenuate each other (as in the decaying periodic function above).

### Multiple time-varying components

A slow-varying component can be added to a rapid component:

`./gaplac sample "y :~| Cat(k) * (SExp(x; l=2) + OU(x; l=0.5) * 1(0.1))" --at "k=1:4;x/k=-5:0.1:5" --plot sample_plot.png --plotx x:k`

### Hierarchical group effects

Suppose you have an experiment with mice. Mice are often grouped into cages, and there is often a very strong cage effect. We can model this as a hierarchical model with one component describing the changes within a cage, and another component describing the individual mouse's differences from the cage mean:

`./gaplac sample "y :~| Cat(cage) * (SExp(x; l=2.5) + Cat(mouse) * SExp(x; l=0.3) * Constant(0.05))" --at "mouse=1:9;cage=floor.((mouse.-1)./3);x/mouse=-5:0.1:5" --plot sample_plot.png --plotx x:mouse`

### "Global" effects

## Variance parameters

***discuss how variance parameters are handled***

## Data likelihoods

### Predictions



