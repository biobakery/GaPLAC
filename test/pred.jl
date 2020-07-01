using DataFrames, CSV, Plots

# julia --project cli/main.jl mcmc "bug :~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Linear(nutrient) + Noise" --data test/testin/input_pair_109.tsv --output test/testin/mcmc_109.tsv --samples 100 --burnin 50 --thin 2
# julia --project cli/main.jl mcmc "bug :~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Noise" --data test/testin/input_pair_109.tsv --output test/testout/mcmc2_109.tsv --samples 100 --burnin 50 --thin 2
# julia --project cli/main.jl select --mcmc test/testin/mcmc_109.tsv --mcmc2 test/testout/mcmc2_109.tsv > test/testout/select_109.txt
# julia --project cli/main.jl predict "bug :~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Linear(nutrient) + Noise" --data test/testin/input_pair_109.tsv --mcmc test/testin/mcmc_109.tsv --at "nutrient=-5:0.1:5;PersonID=0;StoolPairs=0" --output test/testout/prediction_109.tsv 

pred = CSV.File("test/testout/prediction_109.tsv") |> DataFrame
data = CSV.File("test/testin/input_pair_109.tsv") |> DataFrame

x = pred.nutrient
ymu = pred.ymu
y95 = pred.yQ950
y05 = pred.yQ050

scatter(data.bug, data.nutrient)
plot!(x, ymu, ribbon=(y95 .- ymu, ymu .- y05))
xlims!(extrema(data.bug))

describe(data.bug)

# julia --project cli/main.jl mcmc "bug :~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Linear(nutrient) + Noise" --data test/testin/input_pair_3206.tsv --output test/testin/mcmc_3206.tsv --samples 100 --burnin 50 --thin 2
# julia --project cli/main.jl mcmc "bug :~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Noise" --data test/testin/input_pair_3206.tsv --output test/testout/mcmc2_3206.tsv --samples 100 --burnin 50 --thin 2
# julia --project cli/main.jl select --mcmc test/testin/mcmc_3206.tsv --mcmc2 test/testout/mcmc2_3206.tsv > test/testout/select_3206.txt
# julia --project cli/main.jl predict "bug :~| Cat(PersonID) * Cat(StoolPairs) + Cat(PersonID) + Linear(nutrient) + Noise" --data test/testin/input_pair_3206.tsv --mcmc test/testin/mcmc_3206.tsv --at "nutrient=-5:0.1:5;PersonID=0;StoolPairs=0" --output test/testout/prediction_3206.tsv 

pred = CSV.File("test/testout/prediction_3206.tsv") |> DataFrame
data = CSV.File("test/testin/input_pair_3206.tsv") |> DataFrame

x = pred.nutrient
ymu = pred.ymu
y95 = pred.yQ950
y05 = pred.yQ050

scatter(data.bug, data.nutrient)
plot!(x, ymu, ribbon=(y95 .- ymu, ymu .- y05))
xlims!(extrema(data.bug))

