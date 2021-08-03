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

args = Dict()
args["formula"] = "y :~| SExp(:x) * Cat(:subj) + Linear(:diet)"
args["data"] = "gp_sample.csv"
args["output"] = "mcmc_sexp_lin.tsv"


(response, lik, gp_form) = GaPLAC._cli_formula_parse(args["formula"])

df = CSV.read(args["data"], DataFrame)

eq, vars = GaPLAC._apply_vars(gp_form)
kernels = GaPLAC._walk_kernel(eq)
hyper = GaPLAC._hyperparam.(kernels)



GaPLAC._cli_run_mcmc(args)