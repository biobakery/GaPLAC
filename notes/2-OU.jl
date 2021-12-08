using GaPLAC
using AbstractGPs
using KernelFunctions
using Distributions
using CSV
using DataFrames
using Turing
using CairoMakie

(resp, lik, gp_form) = GaPLAC._cli_formula_parse("y ~| OU(:x)")

eq, vars = GaPLAC._apply_vars(gp_form)
kernels = GaPLAC._walk_kernel(eq)


##

se = CSV.read("mcmc.tsv", DataFrame)
ou = CSV.read("mcmc_ou.tsv", DataFrame)
names(se)

fig = Figure()
se_ax = Axis(fig[1,1])
ou_ax = Axis(fig[1,2])