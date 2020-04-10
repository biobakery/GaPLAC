using GPTool
using GPTool.DataFrames
using GPTool.CSV
using Test
using Random

for f in ["gp.pdf", "out.tsv", "sampleplot.png"]
    isfile(f) && rm(f)
end

@testset "BF" begin
end

@testset "Chains" begin
    chains1 = GPTool.Chains()
    @test GPTool.record!(chains1, :a, 1.) == 1
    @test size(chains1.df) == (1,1)
    @test chains1.df.a == [1.]
    
    df = DataFrame(a=rand(10))
    chains2 = GPTool.Chains(df)
    @test GPTool.getrecords(chains2, :a) == df.a

    @test GPTool.record!(chains2, :a, 1.) == 1.
    @test chains2.df.a[10] == 1.
    @test GPTool.record!(chains2, :b, 1.) == 1.
    @test chains2.df.b[10] == 1.
    @test all(isnan, chains2.df.b[1:9])

    @test size(GPTool.thin(chains2, 3, 2).df) == (4,2)

end


@testset "Direct GP" begin
end


@testset "Forumla" begin
end


@testset "GP" begin
end


@testset "Laplace GP" begin
end


@testset "MCMC" begin
end


@testset "MCMC GP" begin
end


@testset "Vis" begin
end


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

    # # need better tests since random generation is not considered stable accross julia versions
    # @test isapprox(first(f), 0.1393494648232693, atol=1e-5)
    # @test isapprox(first(y), 0.14029547194060588, atol=1e-5)

    gp = GPTool.parse_gp_formula("y*Reads/100 : Binomial(Reads) ~| Cat(Person) * SExp(Time) + Noise", ["Time", "y", "Reads", "Person"], [false, false, false, true])
    @test gp.xfun_params == [:Person, :Time]
    @test gp.θc == [0.1,1.,1.]
end


# @testset "Old Tests" begin
#     include("oldtests.jl")
#     for f in ["gp.pdf", "out.tsv", "sampleplot.png"]
#         @test isfile(f)
#         rm(f)
#     end
#
#     dif = run(`diff out.tsv testout/seed1_out.tsv`) |> read |> String
#     @test dif == ""
# end
