module GaPLAC

export invnormaltransform,
       gp_spec,
       likelihood,
       response,
       formula

# GP Packages
using AbstractGPs
using KernelFunctions
using GPLikelihoods
using Distances: Euclidean

# Bayesian stuff
using Turing
using Distributions
using Turing: Variational

# Data stuff
using CSV
using DataFrames
using StatsBase

# Other stuff
using CairoMakie
using ReTest

@testset "Module" begin
    @test true
end

include("gp_parts.jl")
include("liklihoods.jl")
include("abstractgp_translations.jl")
include("utils.jl")
include("interface.jl")
include("plotting.jl")

end # module
