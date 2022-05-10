using CLI
using LoggingExtras

args = CLI.parse_cmdline()

if args["debug"]
    loglevel = Logging.Debug
elseif args["verbose"]
    loglevel = Logging.Info
elseif args["quiet"]
    loglevel = Logging.Error
else
    loglevel = Logging.Warn
end

CLI.setup_logs!(loglevel, args["log"])
@info "Getting started!"

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
