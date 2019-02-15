
include("src/chains.jl")
include("src/gp.jl")
include("src/mcmc.jl")
include("src/laplacegp.jl")
include("src/formula.jl")
include("src/mcmcgp.jl")
include("src/bf.jl")

gp = parse_gp_formula("y : Gaussian(.01) ~| SExp(t; 1) * 1(2)", ["t", "y"], [false, false])

t = collect(range(-5.0, 5.0, length=401))
f, y = samplegp(gp.gp, [], [], [], [], t, t)

using Plots
plot(t, y)

gp = parse_gp_formula("y*Reads/100 : Binomial(Reads) ~| Cat(Person) * SExp(Time) + Noise", ["Time", "y", "Reads", "Person"], [false, false, false, true])
