include("../cli/main.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Debug))

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


gp = parse_gp_formula("y : Gaussian(.01) ~| Cat(person) * SExp(time)", ["person", "time"], [false, false])
x, z, y = gp_inputs(gp, DataFrame(person=1:4, time=5:8))

cmd_sample(Dict(
    "data" => nothing,
    "atdata" => nothing,
    "mcmc" => nothing,
    "output" => "stdout",
    "formula" => "y :~| 1(1) * Cat(person) * SExp(time; l=1)",
    "at" => "person=1:3;time/person=range(-5,5,length=6)",
    "plot" => "sampleplot.png",
    "plotx" => "time:person"))



parsedgp = parse_gp_formula("y : Gaussian(.1) ~| 1", ["y", "person", "time"], [false, false, false])

fw, ∇2lπ_fw, lπ = laplace_approx(parsedgp.gp, LowerTriangular(Matrix{Float64}(I,(2,2))),
    [0.5, -0.5], zeros(2,0), [])


parsedgp = parse_gp_formula("y ~| SExp(x)", ["y", "x"], [false, false])
x, z, y = gp_inputs(parsedgp, DataFrame(x=-1:0.5:1, y=rand.([Uniform() for x in -1:0.5:1])))

chains = mcmcgp(parsedgp.gp, x, y, z, [1., 1., 1.], 20)
