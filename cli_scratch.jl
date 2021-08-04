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
GaPLAC._convert2eq(gp_form; hyperparams=Dict(:x=> 2))

l = 2

df = CSV.read(args["data"], DataFrame)

eq, vars = GaPLAC._apply_vars(gp_form, hyperparams=Dict(:x=>2))
kernels = GaPLAC._walk_kernel(eq)
hyper = GaPLAC._hyperparam.(kernels)


eq

y = df[!, response]
x = Matrix(df[!, vars])

k = GaPLAC.CategoricalKernel() ⊗ SqExponentialKernel()
k2 = with_lengthscale(k, 2)



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

mutable struct Leaf
    name
end

mutable struct Node
    nodes
    properties
end

tree = (:join, 
    (:split,
        "l1",
        (:join, 
            "l2",
            "l3"
        )
    ),
    (:split,
        "l4",
        "l5"
    )
)

props = ["a", "b", "c", "d", "e"]

solution = [
    Node([Node([Leaf("l1")], ["a"]),
          Node([Leaf("l2"), Leaf("l3")], ["b", "c"]),
          ], ["a", "b", "c"]),
    Node([Node([Leaf("l4")], ["d"]),
          Node([Leaf("l5")], ["e"])
          ], ["d", "e"])
]


nested_length(t) = 1
nested_length(t::Tuple) = nested_length(t[2]) + nested_length(t[3])
nested_length(t::Tuple{T}) where T <: Kernel = nested_length.(t.kernels)
nested_length(t::KernelFunctions.TransformedKernel) = nested_length(t.kernel)
nested_length(t::KernelFunctions.KernelTensorProduct) = sum(nested_length.(t.kernels))
nested_length(t::KernelFunctions.KernelSum) = sum(nested_length.(t.kernels))


nested_length.(eq.kernels)
nested_length.(eq.kernels[1])
nested_length(tree)