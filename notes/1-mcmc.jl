using GaPLAC
using AbstractGPs
using KernelFunctions
using Distributions
using CSV
using DataFrames
using Turing

(resp, lik, gp_form) = GaPLAC._cli_formula_parse("y ~| SExp(:x)")

df = CSV.read("data.tsv", DataFrame)
eq, vars = GaPLAC._apply_vars(gp_form)

kernels = GaPLAC._walk_kernel(eq)
inferable = [:x]

y = df[!, resp]
x = Matrix(df[!, vars])

@model function inference_engine(Y, X, eq, inferable)
    ℓ ~ TruncatedNormal(1, 5, 0, 30)
    gp, = GaPLAC._apply_vars(eq; hyperparams=Dict(v=> ℓ for v in inferable))
    
    fx ~ AbstractGPs.FiniteGP(GP(gp), RowVecs(X), 0.1)
    Y .~ Normal.(fx, 1)
end

m = inference_engine(y, x, gp_form, inferable)
    
chain = sample(m, NUTS(0.65), 200, progress=true)
