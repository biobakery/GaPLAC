module GaPLAC

using AbstractGPs
using KernelFunctions
using Distances
using Distributions
using CSV
using DataFrames
using CairoMakie
using Turing
using Turing: Variational
using StatsFuns
using Flux

include("gps.jl")
include("liklihoods.jl")
include("abstractgp_translations.jl")
include("utils.jl")

include("cli/formula.jl")
include("cli/mcmc.jl")
include("cli/formula.jl")

end # module
