using GaPLAC
using AbstractGPs
using KernelFunctions
using Distributions
using CSV
using DataFrames
using Turing
using CairoMakie

spec1 = GaPLAC.gp_spec("y ~| OU(:x)")
spec2 = GaPLAC.gp_spec("y ~| OU(:t) * Categorical(:g) + Linear(:d)")
spec3 = GaPLAC.gp_spec("y ~| OU(:x)")

eq, vars = GaPLAC._apply_vars(GaPLAC.formula(spec))
kernels = GaPLAC._walk_kernel(eq)


##

se = CSV.read("mcmc.tsv", DataFrame)
ou = CSV.read("mcmc_ou.tsv", DataFrame)
names(se)

fig = Figure()
se_ax = Axis(fig[1,1])
ou_ax = Axis(fig[1,2])