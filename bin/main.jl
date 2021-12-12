using ArgParse
using LoggingExtras
using TerminalLoggers

function parse_cmdline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "mcmc"
            help = "Run MCMC to optimize hyperparameters"
            action = :command
        "predict"
            help = "Calculate the posterior of a GP at tdata given data"
            action = :command
        "sample"
            help = "Sample the posterior of a GP"
            action = :command
        "fitplot"
            help = "Diagnostic plots showing the posteriors of different components of the GP"
            action = :command
        "select"
            help = "Output model selection parameters"
            action = :command
        "--verbose", "-v"
            help = "Log level to @info"
            action = :store_true
        "--quiet", "-q"
            help = "Log level to @warning"
            action = :store_true
        "--debug"
            help = "Log level to @debug"
            action = :store_true
        "--log"
            help = "Log to a file as well as stdout"

    end

    @add_arg_table! s["sample"] begin
        "spec"
            help = "GP formula specification"
            required = true
        "--at"
            help = "Range to sample at, eg 'x=-5:0.1:5"
            required = true
        "--plot"
            help = "File to plot to"
        "--output", "-o"
            help = "Table output of GP sample - must end with '.csv' or '.tsv'"
    end

    @add_arg_table! s["mcmc"] begin
        "formula"
            help = "GP formula specification"
            required = true
        "--data", "-i"
            required = true
            help = """
                Table input on which to run inference.
                Must contain columns that correspond to values in 'formula'
                """
        "--infer"
            help = """
                Which model hyperparameter to infer.
                Specify variable names, the hyperparameter(s) will be determined
                based on kernel type (eg length scale for SExp)
                """
            nargs = '+'
            required = true
        "--samples"
            help = """
            Number of samples to take from the chain.
            Default=200
            """
            arg_type=Int
            default=200
        "--output", "-o"
            help = "Table to output sampling chain"
        "--plot"
            help = "File to plot to"
    end

    @add_arg_table! s["select"] begin
        "--formulae"
            help = """
                Compare 2 GP formula specifications, requires '--data' as well.
                Result will be logpdf of formula 2 - logpdf of formula 1.

                A positive value indicates more evidence for formula 2.
                """
            nargs = 2
        "--chains"
            help = """
                Compare 2 sampling chains from 'mcmc' command.
                Result will be the log2 bayes factor.

                A positive value indicates more evidence for chain 1.
                """
            nargs = 2
        "--data", "-i"
            help = """
                Table input on which to run inference.
                Must contain columns that correspond to values in both 'formulae'
                """
        "--plot"
            help = "File to plot to"
    end


    return parse_args(s)
end


# Main
args = parse_cmdline()

function setup_logs!(loglevel, logpath; dryrun=false)
    glog = TerminalLogger(stderr, loglevel)
    if logpath === nothing || dryrun
        global_logger(glog)
    else
        logpath = abspath(expanduser(logpath))
        global_logger(
            TeeLogger(
                MinLevelLogger(FileLogger(logpath), loglevel),
                glog))
    end
end

if args["debug"]
    loglevel = Logging.Debug
elseif args["verbose"]
    loglevel = Logging.Info
elseif args["quiet"]
    loglevel = Logging.Error
else
    loglevel = Logging.Warn
end

setup_logs!(loglevel, args["log"])
@info "Getting started!"

include("cli/GaPLAC_CLI.jl")
using .CLI

args["%COMMAND%"] == "sample" && CLI.Sample.run(args["sample"])
args["%COMMAND%"] == "mcmc" && CLI.MCMC.run(args["mcmc"])
args["%COMMAND%"] == "select" && CLI.Select.run(args["select"])

nothing
# # Redirect logging to a file if necessary
# loglevel = args["debug"] ? Logging.Debug : Logging.Info
# if !isa(args["log"], Nothing)
#     io = open(args["log"], "w+")
#     logger = SimpleLogger(io, loglevel)
#     global_logger(logger)
# else
#     logger = ConsoleLogger(stderr, loglevel; show_limited=false)
#     global_logger(logger)
# end

# # Process commands

# # Clean up
# if !isa(args["log"], Nothing)
#     flush(io)
#     close(io)
# end
