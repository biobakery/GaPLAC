using GPTool
using GPTool.DataFrames
using GPTool.CSV
using Test
using Random

@testset "New stuff"
    @testset "Formula parsing"
        form1 = "y :~| SExp(x; l=1)"
        form2 = "y : Gaussian(.01) ~| SExp(t; l=1.5)*Constant(1)"
        form3 = "bug :~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Linear(nutrient) + Noise"

        (y1, cov1, right1) = GPTool.splitformula(form1)
        (y2, cov2, right2) = GPTool.splitformula(form2)
        (y3, cov3, right3) = GPTool.splitformula(form3)

        
    end
end #new stuff
# @testset "Jason Stuff" begin
#     @testset "BF" begin
#     end

#     @testset "Chains" begin
#         chains1 = GPTool.Chains()
#         @test GPTool.record!(chains1, :a, 1.) == 1
#         @test size(chains1) == (1,1)
#         @test chains1.df.a == [1.]
        
#         df = DataFrame(a=rand(10))
#         chains2 = GPTool.Chains(df)
#         @test GPTool.getrecords(chains2, :a) == df.a

#         @test GPTool.record!(chains2, :a, 1.) == 1.
#         @test chains2.df.a[10] == 1.
#         @test GPTool.record!(chains2, :b, 1.) == 1.
#         @test chains2.df.b[10] == 1.
#         @test all(isnan, chains2.df.b[1:9])

#         @test size(GPTool.thin(chains2, 3, 2)) == (4,2)

#         GPTool.newsample!(chains2)
#         @test size(chains2) == (11,2)

#         # TODO: add tests for read/write chains

#     end


#     @testset "Direct GP" begin
#     end


#     @testset "Forumla" begin
#         @test GPTool.sanitize_var_name("this is a test") == "this_is_a_test"
#         @test GPTool.sanitize_var_name("this.is.a.test") == "this_is_a_test"
#         @test GPTool.sanitize_var_name("σ") == "_"
#         @test GPTool.sanitize_var_name("1234") == "X1234"
#         @test GPTool.sanitize_var_name("_1234") == "_1234"

#         (vn, vs) = GPTool.get_varset(["this is a test", "this.is.a.test", "σ", "1234", "_1234"])
#         @test vn == ["this_is_a_test", "this_is_a_test", "_", "X1234", "_1234"]
#         @test vs == Set([:this_is_a_test, :_, :X1234, :_1234])


        
#         gp = GPTool.parse_gp_formula("y : Gaussian(.01) ~| SExp(t; l=1.5)*Constant(1)", ["t", "y"], [false, false])
#         @test gp isa GPTool.ParsedGPFormula
#         @test length(gp.xnames) == 1
#         @test gp.xnames[1] == "t"
#     end

#     @testset "GP" begin
#     end


#     @testset "Laplace GP" begin
#     end


#     @testset "MCMC" begin
#     end


#     @testset "MCMC GP" begin
#     end


#     @testset "Vis" begin
#     end


#     @testset "API" begin
        
#         isdir("testout/") && rm("testout", recurssive=true, force = true)
#         mkdir("testout/")
#         Random.seed!(42)

#         gp = GPTool.parse_gp_formula("y : Gaussian(.01) ~| SExp(t; l=1.5)*Constant(1)", ["t", "y"], [false, false])
#         @test gp isa GPTool.ParsedGPFormula
#         @test length(gp.xnames) == 1
#         @test gp.xnames[1] == "t"
#         # what else ??
        
#         t = collect(-5.:0.1:5.)
#         f, y = GPTool.samplegp(gp.gp, [], [], [], [], t, t)
#         @test length(f) == length(y)
        
#         # # need better tests since random generation is not considered stable accross julia versions
#         # @test isapprox(first(f), 0.1393494648232693, atol=1e-5)
#         # @test isapprox(first(y), 0.14029547194060588, atol=1e-5)
        
#         gp = GPTool.parse_gp_formula("y*Reads/100 : Binomial(Reads) ~| Cat(Person) * SExp(Time) + Noise", ["Time", "y", "Reads", "Person"], [false, false, false, true])
#         @test gp.xfun_params == [:Person, :Time]
#         @test gp.θc == [0.1,1.,1.]
        
#     end


#     # @testset "Old Tests" begin
#     #     include("oldtests.jl")

#     #     dif = run(`diff out.tsv testout/seed1_out.tsv`) |> read |> String
#     #     @test dif == ""
        
#     #     for f in ["gp.pdf", "out.tsv", "sampleplot.png"]
#     #         @test isfile(f)
#     #         rm(f)
#     #     end
#     # end

#     @testset "CLI" begin
#         pred = "testout/prediction.tsv"
#         isfile(pred) && rm(pred)
#         run(`julia ../cli/main.jl predict "bug ~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Linear(nutrient)" --data ./testin/input_pair_1003.tsv --mcmc ./testin/mcmc_1003.tsv --at "nutrient=-5:0.1:5;PersonID=0;StoolPairs=0" --output testout/prediction.tsv`)
#         @test isfile(pred)
#     end
# end #jason stuff