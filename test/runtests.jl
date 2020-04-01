using GPTool
using GPTool.DataFrames
using GPTool.CSV
using Test

@testset "API" begin
    
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
