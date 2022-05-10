module CLI

using ArgParse
using LoggingExtras
using TerminalLoggers
using Statistics

include("mcmc.jl")
include("sample.jl")
include("select.jl")
include("main.jl")

end # module