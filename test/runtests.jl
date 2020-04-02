using GPTool
using GPTool.DataFrames
using GPTool.CSV
using Test
using Random

@testset "API" begin
    Random.seed!(42)

    gp = GPTool.parse_gp_formula("y : Gaussian(.01) ~| SExp(t; l=1.5)*Constant(1)", ["t", "y"], [false, false])
    @test gp isa GPTool.ParsedGPFormula
    @test length(gp.xnames) == 1
    @test gp.xnames[1] == "t"
    # what else ??

    t = collect(-5.:0.1:5.)
    f, y = GPTool.samplegp(gp.gp, [], [], [], [], t, t)
    @test length(f) == length(y)
    @test first(f) == 0.1393494648232693
    @test first(y) == 0.14029547194060588

    gp = GPTool.parse_gp_formula("y*Reads/100 : Binomial(Reads) ~| Cat(Person) * SExp(Time) + Noise", ["Time", "y", "Reads", "Person"], [false, false, false, true])
    @test gp.xfun_params == [:Person, :Time]
    @test gp.θc == [0.1,1.,1.]
end


@testset "Old Tests" begin
    @info pwd()
    for f in ["gp.pdf", "out.tsv", "sampleplot.png"]
        isfile(f) && rm(f)
    end

    include("oldtests.jl")

    dif = run(`diff out.tsv testout/seed1_out.tsv`) |> read |> String
    @test dif == ""
end
