using CSV
using DataFrames
using Turing
using Distributions
using StatsPlots
using Distances
using LinearAlgebra
using StatsModels
using AbstractGPs, KernelFunctions


df = DataFrame(
    subject = repeat(1:5, inner=3),
    obese = repeat(rand(Bool, 5), inner=3),
    timepoint = [1,2,3,1,3,4,1,2,5,1,4,5,1,3,5],
    bug = rand(Beta(0.9, 5), 15),
    nutrient = rand(Beta(0.9,5), 15)
)

ip1 = CSV.File("test/testin/input_pair_109.tsv") |> DataFrame

@model GPRegression(y, df) = begin
    param_1, param_2, param_3, noise = ones(4)
    k_nutrient = AbstractGPs.transform(LinearKernel(), param_1) # kernel for nutrient
    k_obese = AbstractGPs.transform(LinearKernel(), param_2) # kernel for obese
    k_t = AbstractGPs.transform(SqExponentialKernel(), param_3) # kernel for time
    k = TensorProduct([k_nutrient, k_obese, k_t]) # Collect all the kernels to make them act dimension wise
    
    # Here we create a the prior based on the kernel and the data
    priorgp = AbstractGPs.FiniteGP(GP(k), hcat(df.obese, df.nutrient, df.timepoint), noise, obsdim = 1) # Not sure about the dimensionality here
    f ~ priorgp
    for i in 1:size(df,1)
        y[i] ~ Normal(f[i], noise)
    end

end

gp = GPRegression(df.bug, df)

@time chain = sample(gp, HMC(0.01, 100), 200)

using Turing, Distributions

# Import MCMCChains, Plots, and StatPlots for visualizations and diagnostics.
using MCMCChains, StatsPlots

# Functionality for splitting and normalizing the data.
using MLDataUtils: shuffleobs, splitobs, rescale!

# Functionality for evaluating the model predictions.
using Distances