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

sqexpkernel(alpha::Real, rho::Real) = 
    alpha^2 * AbstractGPs.transform(SqExponentialKernel(), 1/(rho*sqrt(2)))

# Define model.
@model GPRegression(y, subj, x) = begin
    nsub = length(unique(subj))
    # GP Priors.
    alpha ~ LogNormal(0.0, 0.1)
    rho ~ LogNormal(0.0, 1.0)
    gp_sigma ~ LogNormal(0.0, 1.0)
    
    lin_sigma ~ LogNormal(0.0, 1.0)
    # Covariance function.
    kernel = sqexpkernel(alpha, rho)
    # GP (implicit zero-mean)
    gp = GP(kernel)
    
    # index variable for subject
    si ~ filldist(Normal(), length(unique(subj)))
    mu = si[subj]
    
    # Sampling Distribution (MvNormal likelihood).
    y .~ Normal.(mu, lin_sigma)
end;


m = GPRegression(ip1.bug, ip1.PersonID, ip1.Date, ip1.nutrient)

# Fit via HMC.
@time chain = sample(m, HMC(0.01, 100), 200)  # start sampling.

plot(chain)


x = vcat(randn(10), randn(10).+10, randn(10).-5)
y = vcat(randn(10), randn(10).+10, randn(10).-5) .* 2
g = vcat(ones(Int, 10), ones(Int, 10)*2, ones(Int, 10)*3)
tp = repeat(1:10, 3)

sekernel(alpha, rho) = 
  alpha^2 * KernelFunctions.transform(SEKernel(), sqrt(0.5)/rho)

@model function demo(x, y, g, tp, jitter=1e-6)
    ng = length(unique(g))
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

model = demo(x,y,g,tp);
sample(model, NUTS(0.8), 1_000)

pidmap = Dict(p=>i for (i,p) in enumerate(unique(ip1.PersonID)))
ip1.pid = [pidmap[p] for p in ip1.PersonID]


md2 = demo(ip1.nutrient[1:50], ip1.bug[1:50], ip1.pid[1:50], ip1.Date[1:50])
r = sample(md2, NUTS(0.8), 1_000)
md3 = demo(randn(50), ip1.bug[1:50], ip1.pid[1:50], ip1.Date[1:50])
r2 = sample(md3, NUTS(0.8), 1_000)

