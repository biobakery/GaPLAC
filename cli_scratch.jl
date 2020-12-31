include("dietmodels.jl")

infile = ARGS[1]
outdir = ARGS[2]
prefix = ARGS[3]

@info "Getting started" infile outdir prefix

!isdir(outdir) && mkdir(outdir)

df = CSV.read(infile, DataFrame)
df = df[completecases(df), :]
pidmap = Dict(p=>i for (i,p) in enumerate(unique(df.PersonID)))
df.pid = [pidmap[p] for p in df.PersonID]

df.datemod = [d+rand(Normal(0, 0.2)) for d in df.Date]
df.bugmod = [b+rand(Normal(0, 0.2)) for b in df.bug]

scorr = subjectcorrmat(df.pid)

@info "Sampling from first model"

gpm1 = GPmodel1(df.bugmod, df.nutrient, scorr, df.datemod)
r1 = sample(gpm1, HMC(0.1,20), 100)

@info "Sampling from second model"

gpm2 = GPmodel2(df.bugmod, scorr, df.datemod)
r2 = sample(gpm2, HMC(0.1,20), 100);

@info "Writing outputs"

l2b = log2bayes(r1, r2)

baysfile = joinpath(outdir, prefix*"_log2bayes.txt")
write(bayesfile, string(l2b)*'\n')

result_df = DataFrame()
result_df.model1_lp = r1[:lp][chain=1] |> collect
result_df.model1_σ2 = r1[:σ2][chain=1] |> collect
result_df.model1_c = r1[:c][chain=1] |> collect 
result_df.model2_lp = r2[:lp][chain=1] |> collect
result_df.model2_σ2 = r2[:σ2][chain=1] |> collect

result_file = joinpath(outdir, prefix*"_results.csv")
CSV.write(result_file, result_df)