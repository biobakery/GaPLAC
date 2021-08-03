using GaPLAC
using KernelFunctions
using CairoMakie
using AbstractGPs

args = Dict()
args["formula"] = "y :~| SExp(:x; l=2.5) * Cat(:subj) + Linear(:diet)"
args["at"] = "x=-4:0.5:4; diet=-2:1:2"
args["output"] = "gp_sample.csv"
args["plot"] = "gp_sample.png"

df, gp = GaPLAC._cli_run_sample(args)

## Turing

using GaPLAC

args = Dict()
args["formula"] = "y :~| SExp(:x)"
args["data"] = "gp_sample.csv"
args["output"] = "mcmc_sexp.tsv"

GaPLAC._cli_run_mcmc(args)

##

using GaPLAC
using KernelFunctions
using CSV
using DataFrames


args = Dict()
args["formula"] = "y :~| SExp(:x) * Cat(:subj) + Linear(:diet)"
args["data"] = "gp_sample.csv"
args["output"] = "mcmc_sexp_lin.tsv"

(response, lik, gp_form) = GaPLAC._cli_formula_parse(args["formula"])

gp_form

df = CSV.read(args["data"], DataFrame)

eq, vars = GaPLAC._apply_vars(gp_form)
kernels = GaPLAC._walk_kernel(eq)
hyper = GaPLAC._hyperparam.(kernels)

y = df[!, response]
x = Matrix(df[!, vars])

k = GaPLAC.CategoricalKernel() * SqExponentialKernel()
k2 = with_lengthscale(k, 2)
k3 = GaPLAC.CategoricalKernel() * with_lengthscale(SqExponentialKernel(), 2)

# GaPLAC._cli_run_mcmc(args)

##

using Distributions
using CairoMakie
using StatsFuns

function invnormaltransform(v; μ=0, σ=1, c=3/8)
    rank = invperm(sortperm(v))
    N = length(v)
    return [norminvcdf(μ, σ, (x - c) / (N - 2c + 1)) for x in rank]
end

x = rand(Beta(1,3), 1000); y = invnormaltransform(x);

hist(x)
hist!(y)
current_figure()
