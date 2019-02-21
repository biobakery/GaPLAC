
include("src/chains.jl")
include("src/gp.jl")
include("src/mcmc.jl")
include("src/laplacegp.jl")
include("src/formula.jl")
include("src/mcmcgp.jl")
include("src/bf.jl")

gp = parse_gp_formula("y : Gaussian(.01) ~| SExp(t; l=1.5)*Constant(1)", ["t", "y"], [false, false])

t = collect(-5.:0.1:5.)
f, y = samplegp(gp.gp, [], [], [], [], t, t)

using DataFrames
using CSV
CSV.write("out.tsv", DataFrame(t=t, y=y), delim='\t')

using Plots
plot(t, y, size = (300, 200))
savefig("./gp.pdf")

gp = parse_gp_formula("y*Reads/100 : Binomial(Reads) ~| Cat(Person) * SExp(Time) + Noise", ["Time", "y", "Reads", "Person"], [false, false, false, true])


function foo(f)
    f() + 1
end

function watest()
    f = eval(:(()->1))
    return Base.invokelatest(foo, f)
end

watest()
