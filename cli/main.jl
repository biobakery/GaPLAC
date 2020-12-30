using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using GPTool
using ArgParse


function parse_cmdline()
    s = ArgParseSettings()

    @add_arg_table s begin
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
        "--log"
            help = "Direct logging output to a file"
        "--debug"
            help = "Enable debug logging"
            action = :store_true
    end

    @add_arg_table s["mcmc"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;#subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            required = true
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "none"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters; if provided, the chain will be extended. Does NOT support --data format flags"
        "--burnin", "-r"
            help = "Burn-in samples"
            arg_type = Int
            default = 0
        "--thin", "-t"
            help = "Thinning"
            arg_type = Int
            default = 1
        "--samples", "-n"
            help = "Number of MCMC samples to generate (after anneal, burn-in, and thinning)"
            arg_type = Int
            default = 100
        "--output", "-o"
            help = "Output filename, accepts \"stdout\", does NOT support --data format flags"
            default = "stdout"
    end

    @add_arg_table s["predict"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;#subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            required = true
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "none"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters, does NOT support --data format flags"
        "--atdata", "-t"
            help = "Data files providing variable values at which to make predictions"
        "--at"
            help = "Predict at these points (alternative to tdata); format: ;-separated variable=start:step:end OR variable=value OR variable=Julia distribution OR variable/group=value"
        "--output", "-o"
            help = "Output filename, accepts \"stdout\", supports --data format flags"
            default = "stdout"
    end

    @add_arg_table s["sample"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;#subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "none"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters, does NOT support --data format flags"
        "--atdata", "-t"
            help = "Data files providing variable values at which to sample the GP"
        "--at"
            help = "Predict at these points (alternative to tdata); format: ;-separated variable=start:step:end OR variable=value OR variable=Julia distribution OR variable/group=value"
        "--output", "-o"
            help = "Output filename, accepts \"stdout\", supports --data format flags"
            default = "stdout"
        "--plot"
            help = "Output filename for the sample plots"
        "--plotx"
            help = "X axis of the sample plot; can also be \"x:group\" to provide a grouping variable as well"
    end

    @add_arg_table s["fitplot"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;#subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            required = true
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "none"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters, does NOT support --data format flags"
        "--component", "--comp"
            help = "Plot fit plots for components with the given variables"
        "--output", "-o"
            help = "Output filename"
            default = "fitplots.pdf"
    end

    @add_arg_table s["select"] begin
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters, does NOT support --data format flags"
        "--mcmc2", "-c"
            help = "MCMC samples for hyperparameters of the second model (for model comparison), does NOT support --data format flags"
    end

    return parse_args(s)
end

include("cli_include.jl")

function main(args)
    # Redirect logging to a file if necessary
    loglevel = args["debug"] ? Logging.Debug : Logging.Info
    if !isa(args["log"], Nothing)
        io = open(args["log"], "w+")
        logger = SimpleLogger(io, loglevel)
        global_logger(logger)
    else
        logger = ConsoleLogger(stderr, loglevel; show_limited=false)
        global_logger(logger)
    end

    # Process commands
    processcmd(args)

    # Clean up
    if !isa(args["log"], Nothing)
        flush(io)
        close(io)
    end
end

args = parse_cmdline()

main(args)
