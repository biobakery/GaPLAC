include("../cli/cli_include.jl")

using GaPLAC
using Logging
using Random

Random.seed!(1)

global_logger(ConsoleLogger(stderr, Logging.Info))

gp = GaPLAC.parse_gp_formula("y : Gaussian(.01) ~| SqExp(t; l=1.5)*Constant(1)", ["t", "y"], [false, false])


t = collect(-5.:0.1:5.)
f, y = GaPLAC.samplegp(gp.gp, [], [], [], [], t, t)

using DataFrames
using CSV
CSV.write("out.tsv", DataFrame(t=t, y=y), delim='\t')

using Plots
plot(t, y, size = (300, 200))
savefig("./gp.pdf")

gp = GaPLAC.parse_gp_formula("y*Reads/100 : Binomial(Reads) ~| Cat(Person) * SqExp(Time) + Noise", ["Time", "y", "Reads", "Person"], [false, false, false, true])


function foo(f)
    f() + 1
end

function watest()
    f = eval(:(()->1))
    return Base.invokelatest(foo, f)
end

watest()


gp = GaPLAC.parse_gp_formula("y : Gaussian(.01) ~| Cat(person) * SqExp(time)", ["person", "time"], [false, false])
x, z, y = GaPLAC.gp_inputs(gp, DataFrame(person=1:4, time=5:8))

cmd_sample(Dict(
    "data" => nothing,
    "atdata" => nothing,
    "mcmc" => nothing,
    "output" => "stdout",
    "formula" => "y :~| 1(1) * Cat(person) * SqExp(time; l=1)",
    "at" => "person=1:3;time/person=range(-5,5,length=6)",
    "plot" => "sampleplot.png",
    "plotx" => "time:person"))



parsedgp = GaPLAC.parse_gp_formula("y : Gaussian(.1) ~| 1", ["y", "person", "time"], [false, false, false])

fw, ∇2lπ_fw, lπ = GaPLAC.laplace_approx(parsedgp.gp, GaPLAC.LowerTriangular(Matrix{Float64}(GaPLAC.I,(2,2))),
    [0.5, -0.5], zeros(2,0), [])


parsedgp = GaPLAC.parse_gp_formula("y ~| SqExp(x)", ["y", "x"], [false, false])
x, z, y = GaPLAC.gp_inputs(parsedgp, DataFrame(x=-1:0.5:1, y=rand.([Uniform() for x in -1:0.5:1])))

chains = GaPLAC.mcmcgp(parsedgp.gp, x, y, z, [1., 1., 1.], 20)
