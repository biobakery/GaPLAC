module GPTool

using DataFrames



include("src/chains.jl")
include("src/gp.jl")
include("src/mcmc.jl")
include("src/laplacegp.jl")
include("src/directgp.jl")
include("src/formula.jl")
include("src/mcmcgp.jl")
include("src/bf.jl")

end # module
