using CSV
using DataFrames
using Turing
using Distributions
using StatsPlots
using Distances
using LinearAlgebra
using AbstractGPs, KernelFunctions


df = DataFrame(
    subject = repeat(1:5, inner=3),
    obese = repeat(rand(Bool, 5), inner=3),
    timepoint = [1,2,3,1,3,4,1,2,5,1,4,5,1,3,5],
    bug = rand(Beta(0.9, 5), 15),
    nutrient = rand(Beta(0.9,5), 15)
)

ip1 = CSV.File("test/testin/input_pair_109.tsv") |> DataFrame

function numbsubj(subjects)
    id = Dict(s=> i for (i,s) in enumerate(subjects))

    return [id[s] * 1000 for s in subjects]
end

ip1.subject = numbsubj(ip1.PersonID)


sekernel(alpha, rho) = 
  alpha^2 * KernelFunctions.transform(SEKernel(), sqrt(0.5)/rho)

@model function myGP(y, X, T, jitter=1e-6)
    N, P = size(T)  # Note that T needs to have dimension Nx1 (a matrix)
    
    # Dimensions of linear model predictors
    J = size(X, 2)  # X should be N x J
    
    # Priors.
    mu ~ Normal(0, 1)
    sig2 ~ LogNormal(0, 1)
    alpha ~ LogNormal(0, 0.1)
    rho ~ LogNormal(0, 1)

    beta ~ filldist(Normal(0, 1), J) # Prior for linear model coefficients
    
    # GP Covariance matrix
    kernel = sekernel(alpha, rho)  # covariance function
    K = kernelmatrix(kernel, T, obsdim=1)  # cov matrix
    K += LinearAlgebra.I * (sig2 + jitter) # jitter for numerical stability

    y ~ MvNormal(mu .+ X * beta, K) # Sampling Distribution.
end

gp = myGP(df.bug,
          Matrix(df[!,[:subject, :obese, :nutrient]]),
          Matrix{Float64}(df[!,[:timepoint]]))

@time chain = sample(gp, HMC(0.01, 100), 200)

plot(chain)
savefig("~/Desktop/plot.png")
ip1.t = map(t-> t + randn()/100, ip1.StoolPairs)
ip1.s = map(t-> t + randn()/100, ip1.subject)


gp2 = myGP(ip1.bug,
    Matrix(ip1[!,[:s, :nutrient]]),
    Matrix{Float64}(ip1[!,[:t]]))

@time chain2 = sample(gp2, HMC(0.01, 100), 200)

plot(chain2)
