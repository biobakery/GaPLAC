using CSV
using DataFrames
using Turing
using Distributions
using StatsPlots
using Distances
using AbstractGPs, KernelFunctions


ip1 = CSV.File("test/testin/input_pair_109.tsv") |> DataFrame
pidmap = Dict(p=>i for (i,p) in enumerate(unique(ip1.PersonID)))
ip1.pid = [pidmap[p] for p in ip1.PersonID]


sekernel(alpha, rho) = 
  alpha^2 * KernelFunctions.transform(SEKernel(), sqrt(0.5)/rho)

@model function GPmodel1(x, y, g, tp, jitter=1e-6)
    ng = length(unique(g))
    x = x .+ rand(Normal(0.0, 0.01), length(x))
    y = y .+ rand(Normal(0.0, 0.01), length(y))
    # group coefficients priors
    a ~ filldist(Normal(), ng)
    
    # x coefficient prior
    b ~ Normal()

    # gp priors
    sig2 ~ LogNormal(0, 1)
    alpha ~ LogNormal(0, 0.1)
    rho ~ LogNormal(0, 1)

    # lm specs
    mu = a[g] + b * x
    
    # gp specs
    kernel = sekernel(alpha, rho)  # covariance function
    K = kernelmatrix(kernel, reshape(tp, length(tp), 1), obsdim=1)  # cov matrix
    K += LinearAlgebra.I * (sig2 + jitter)
    
    y .~ MvNormal(mu, K)
end

md1 = GPmodel1(ip1.nutrient[1:100], ip1.bug[1:100], ip1.pid[1:100], ip1.Date[1:100])
r1 = sample(md1, NUTS(0.8), MCMCThreads(), 6, 100)



@model function GPmodel2(y, g, tp, jitter=1e-6)
    ng = length(unique(g))
    y = y .+ rand(Normal(0.0, 0.01), length(y))
    # group coefficients priors
    a ~ filldist(Normal(), ng)
    
    # x coefficient prior
    b ~ Normal()

    # gp priors
    sig2 ~ LogNormal(0, 1)
    alpha ~ LogNormal(0, 0.1)
    rho ~ LogNormal(0, 1)

    # lm specs
    mu = a[g]
    
    # gp specs
    kernel = sekernel(alpha, rho)  # covariance function
    K = kernelmatrix(kernel, reshape(tp, length(tp), 1), obsdim=1)  # cov matrix
    K += LinearAlgebra.I * (sig2 + jitter)
    
    y .~ MvNormal(mu, K)
end

md2 = GPmodel2(ip1.bug[1:100], ip1.pid[1:100], ip1.Date[1:100])
r2 = sample(md2, NUTS(0.8), MCMCThreads(), 6, 100)

