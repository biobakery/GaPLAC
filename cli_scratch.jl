using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include("dietmodels.jl")

infile = ARGS[1]
outdir = ARGS[2]
prefix = ARGS[3]
chain_length = parse(Int, ARGS[4])
# pairid = 3206
# infile = "test/testin/input_pair_$pairid.tsv"
# outdir = "test/testout/"
# prefix = "pair_$pairid-2"

@info "Getting started" infile outdir prefix

!isdir(outdir) && mkdir(outdir)

df = CSV.read(infile, DataFrame)
df = df[completecases(df), :]
pidmap = Dict(p=>i for (i,p) in enumerate(unique(df.PersonID)))
df.pid = [pidmap[p] for p in df.PersonID]
sort!(df, :pid)
dfg = groupby(df, :pid)
grouped_bugs = [grp.bug for grp in dfg]
grouped_tps = [grp.Date for grp in dfg]
grouped_diet = [grp.nutrient for grp in dfg]
pids = unique(df.pid)

@info "Sampling from first model"

gpm1 = GPmodel1(grouped_bugs, grouped_tps, pids)
@time r1 = sample(gpm1, HMC(0.1,20), chain_length)

@info "Sampling from second model"

gpm2 = GPmodel2(grouped_bugs, grouped_diet, grouped_tps, pids)
@time r2 = sample(gpm2, HMC(0.1,20), chain_length)

@info "Writing outputs"

open(joinpath(outdir, prefix*"_log2bayes.txt"), "w") do io
    l2b, lp2, lp1 = log2bayes(r2, r1)
    println(io, "Log2 bayes: ", l2b)
    println(io, "Log(p) without diet: ", lp1)
    println(io, "Log(p) with diet: ", lp2)
end

result_df = DataFrame()
result_df.model1_lp = r1[:lp][chain=1] |> collect
result_df.model1_σ2 = r1[:sigma][chain=1] |> collect

result_df.model2_lp = r2[:lp][chain=1] |> collect
result_df.model2_σ2 = r2[:sigma][chain=1] |> collect
result_df.model2_c = r2[:c][chain=1] |> collect
result_df.model2_m = r2[:m][chain=1] |> collect 

CSV.write(joinpath(outdir, prefix*"_results.csv"), result_df)
CSV.write(joinpath(outdir, prefix*"_model1.csv"), DataFrame(r1)) # no diet
CSV.write(joinpath(outdir, prefix*"_model2.csv"), DataFrame(r2)) # with diet