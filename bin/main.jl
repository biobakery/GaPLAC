using GaPLAC
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
            help = "Output model selection parameters; requires --mcmc and --mcmc2"
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
        "formula"
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
                Must contain columns that correspond to values in 'formula'x
                """
        "--output", "-o"
            help = "Table to output sampling chaing"
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


args["%COMMAND%"] == "sample" && GaPLAC._cli_run_sample(args["sample"])
args["%COMMAND%"] == "mcmc" && GaPLAC._cli_run_mcmc(args["mcmc"])

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
