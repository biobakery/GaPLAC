
# Guide

[![Build Status](https://github.com/biobakery/gptool.jl/workflows/CI/badge.svg)](https://github.com/biobakery/gptool.jl/actions)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kescobo.github.io/gptool.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kescobo.github.io/gptool.jl/dev)

This guide is intended to provide an overview of the basic workflow using GaPLAC. A complete command reference, as well as available covariance and likelihood functions are provided below.

## Installation

1. Install [Julia](https://julialang.org/).
2. Download GaPLAC's repository and unpack it somewhere.
3. Open a console in GaPLAC's root folder and run 
   ```
   $ julia --project=@. -e 'using Pkg; Pkg.instantiate()'`

   This will install the required packages and may take a few minutes.
4. If on Mac/UNIX, to use `./gaplac ...` format, you may need to run `chmod u+x ./gaplac`.

## Generating some sample data

GaPLAC has five main commands to work with: `sample`, `mcmc`, `select`, `predict`, and `fitplot`. We will first look at `sample`, which draws a sample from a Gaussian Process. This can be helpful to visualize the kinds of functions described by a particular GP, or to provide some sample data to test later functions.

Run the following command from the GaPLAC root folder:

```
./gaplac sample "y :~| SqExp(:x; l=1)" --at "x=-5:0.1:5" --plot gp_sample.png
```

This may take a few minutes the first time, since Julia must compile all the packages. It will produce a large amount of output to the console, and should also produce a plot in `gp_sample.png` which resembles:

![Wavy line](img/guide1.png)

If instead you get an error mentioning missing dependencies, it means that step 4 of the Installation section was not successfully completed. Try following the Installation instructions again to resolve this.

Let's look at each of the pieces of the command:

- `"y :~| SqExp(:x; l=1)"`: This is the GP formula, much like a model formula in R. In this case, the output (`y`) is modeled as a GP with a Squared-Exponential covariance function (`SqExp`) with a lengthscale (`l`) of `1`. Note also the `:` in `:~|`. Normally, a data likelihood (described later) can be specified between the `:` and the `~`, but here we don't specify anything, and the GP will be modeled _without_ a likelihood. This effectively allows us more directly observe the types of dynamics modeled by the Gaussian Process described in the formula.
- `--at "x=-5:0.1:5"`: This tells GaPLAC what values of `x` to sample the GP at.
- `--plot gp_sample.png`: Plot the dynamics here.

Try changing the lengthscale of the `SqExp` term. How does this affect the function? Try adding other components (the full list is at the end of this document) by adding them to the formula, such as an Ornstein-Uhlenbeck process (`OU(x; l=1)`), or simply some `Noise`.

Now let's generate a smaller set of data at some randomly chosen `x` coordinates, and store the results in a file instead of printing to stdout:

```
./gaplac sample "y :~| SqExp(:x; l=1.5)" --at "x = rand(Uniform(-5,5), 50)" --output data.tsv
```

Look at the contents of `data.tsv`. It should contain two columns: `x` and `y`, and the rows are not sorted in any way. We will use this data for the next command.


## Fitting parameters

We are usually interested in the parameters of the covariance function which best fits some data. This is accomplished with the `mcmc` command in GaPLAC, which will produce a [MCMC chain](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) of samples from the posterior distribution of the model parameters.

Try running the following command:

```
./gaplac mcmc "y ~| SqExp(:x)" --data data.tsv --output mcmc.tsv --samples 500 --infer x
```

First, note that the model formula omits the additional `:` before the `~`. This will therefore default to a Gaussian likelihood.

Now examine the output file `mcmc.tsv`. Instead of showing a relationship between `x` and `y`, as the previous run of `sample` did, this file should contain several columns of parameter values, as well as the all-important final column containing the log of the unnormalized posterior density for the sample (think of this something like the goodness of fit). In particular, take a look at the covariance function parameter `ℓ`, which we set above in the model formula to `1.5`, but which we did not tell the `mcmc` command. If all worked well, the mean of this parameter should converge to, and hover around the true value of `1.5`.

## TODO: Add `predict`

## Comparing models

`mcmc` fits a single model, but how do we know it's the right model? Maybe another model would be a better fit? To answer this question, we can run `mcmc` again with the other model. In this case, let's test a different kind of time-varying process called an Ornstein-Uhlenbeck (OU) process:

```
./gaplac mcmc "y ~| OU(:x)" --data data.tsv --output mcmc_ou.tsv --samples 500
```

This will give us a second set of model fit results in a new `mcmc_ou.tsv` file. Now we can use those goodness-of-fit values in each of the files to determine which of the models we believe (hint: we generated data from a Squared-Exponential covariance function, and thus we expect the OU process will perform worse. Let's test that with the `select` command:

```
./gaplac select --chains mcmc.tsv mcmc_ou.tsv
```
```
┌ Info: Log2 Bayes: 8.405
│
│   •  Log(pdf) - model 1: -81.29118
│     
│
│   •  Log(pdf) - model 2: -89.69639
│      9.972996e-28
│
└ Note - Positive values indicate more evidence for model 1
```

This will compare the log posterior values stored in each of the MCMC chains, and summarize them as a [Bayes Factor](https://en.wikipedia.org/wiki/Bayes_factor), which is reported in log2 scale. Here, log2 Bayes Factors greater than 1 indicate that the first model (in this case the Squared-Exponential) should be preferred, while negative numbers indicate the opposite - that the second model should be preferred.

You may also compare different forumla parameters on your initial data,
rather than using the outputs from MCMC.

```
./gaplac -v select --formulae "y ~| SqExp(:x, l=2)" "y ~| SqExp(:x, l=1)" --data data.tsv
```
```
[ Info: running 'select'
┌ Info:
│ Dict{String, Any} with 4 entries:
│   "plot" => nothing
│   "formulae" => Any["y ~| SqExp(:x, l=1.5)", "y ~| OU(:x, l=1.5)"]
│   "data" => "data.tsv"
└   "chains" => Any[]
┌ Info: Log2 Bayes: 4.44
│
│   •  Log(pdf) - model 1: -31.53397005887427
│
│   •  Log(pdf) - model 2: -35.97395926954643
│
└ Note - Positive values indicate more evidence for model 1
```


# Command references

Available by running the scrupt with `--help`

## Commands

```sh
./gaplac --help
usage: main.jl [-v] [-q] [--debug] [--log LOG] [-h]
               {mcmc|predict|sample|fitplot|select}

commands:
  mcmc           Run MCMC to optimize hyperparameters
  predict        Calculate the posterior of a GP given data # not yet implemented
  sample         Sample the posterior of a GP
  fitplot        Diagnostic plots showing the posteriors of different
                 components of the GP # not yet implemented
  select         Output model selection parameters; requires --mcmc
                 and --mcmc2

optional arguments:
  -v, --verbose  Log level to @info
  -q, --quiet    Log level to @warning
  --debug        Log level to @debug
  --log LOG      Log to a file as well as stdout
  -h, --help     show this help message and exit
```

## Sample


```sh
./gaplac sample --help
usage: main.jl sample --at AT [--plot PLOT] [-o OUTPUT] [-h] formula

positional arguments:
  formula              GP formula specification

optional arguments:
  --at AT              Range to sample at, eg 'x=-5:0.1:5
  --plot PLOT          File to plot to
  -o, --output OUTPUT  Table output of GP sample - must end with
                       '.csv' or '.tsv'
  -h, --help           show this help message and exit
```

## MCMC

```sh
./gaplac mcmc --help
usage: main.jl mcmc -i DATA --infer INFER [INFER...] [-o OUTPUT]
                    [--plot PLOT] [-h] formula

positional arguments:
  formula               GP formula specification

optional arguments:
  -i, --data DATA       Table input on which to run inference. Must
                        contain columns that correspond to values in
                        'formula'
  --infer INFER [INFER...]
                        Which model hyperparameter to infer. Specify
                        variable names, the hyperparameter(s) will be
                        determined based on kernel type (eg length
                        scale for SqExp)
  -o, --output OUTPUT   Table to output sampling chain
  --plot PLOT           File to plot to
  -h, --help            show this help message and exit
```

## Select

```sh
./gaplac select --help
usage: main.jl select [--formulae FORMULAE FORMULAE]
                      [--chains CHAINS CHAINS] [-i DATA] [--plot PLOT]
                      [-h]

optional arguments:
  --formulae FORMULAE FORMULAE
                        Compare 2 GP formula specifications, requires
                        '--data' as well. Result will be logpdf of
                        formula 2 - logpdf of formula 1. A positive
                        value indicates more evidence for formula 2.
  --chains CHAINS CHAINS
                        Compare 2 sampling chains from 'mcmc' command.
                        Result will be the log2 bayes factor. A
                        positive value indicates more evidence for
                        chain 1.
  -i, --data DATA       Table input on which to run inference. Must
                        contain columns that correspond to values in
                        both 'formulae'
  --plot PLOT           File to plot to
  -h, --help            show this help message and exit
```
