using GaPLAC
using KernelFunctions
using CairoMakie
using AbstractGPs

args = Dict()
args["formula"] = "y :~| SExp(:x; l=2.5) + Linear(:diet)"
args["at"] = "x=-4:0.5:4"
args["output"] = "gp_sample.csv"
args["plot"] = "gp_sample.png"

df, gp = GaPLAC._cli_run_sample(args)

## Turing

using GaPLAC

args = Dict()
args["formula"] = "y :~| SExp(:x) + Linear(:diet)"
args["data"] = "gp_sample.csv"
args["output"] = "mcmc_sexp.tsv"
args["infer"] = ["x"]

GaPLAC._cli_run_mcmc(args)

##

using GaPLAC
using KernelFunctions
using CSV
using DataFrames
using Turing

args = Dict()
args["formula"] = "y :~| SExp(:x) + Linear(:diet)"
args["data"] = "gp_sample.csv"
args["output"] = "mcmc_sexp_lin.tsv"
args["infer"] = ["x"]

(response, lik, gp_form) = GaPLAC._cli_formula_parse(args["formula"])
@debug "GP formula" gp_form
    
df = CSV.read(args["data"], DataFrame)

eq, vars = GaPLAC._apply_vars(gp_form)
kernels = GaPLAC._walk_kernel(eq)
inferable = Symbol.(args["infer"])

y = df[!, response]
x = Matrix(df[!, vars])

@model function inference_engine(Y, X, eq, inferable)
    ℓ ~ TruncatedNormal(1, 5, 0, 30)
    gp, = GaPLAC._apply_vars(eq; hyperparams=Dict(v=> ℓ for v in inferable))

    fx ~ AbstractGPs.FiniteGP(GP(gp), RowVecs(X), 1e-6)
    Y .~ Normal.(fx, 1)
end

m = inference_engine(y, x, gp_form, inferable)
chain = sample(m, NUTS(0.65), 200, progress=true)

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


##
