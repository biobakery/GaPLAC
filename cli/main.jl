
using CSV

@add_arg_table s begin
    # All commands:
    "formula"
        help = "Gaussian Process formula"
        required = true
    "--data", "-d"
        help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;~subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
        default = "stdin"
    "--bind", "-b"
        help = "Name bindings, format is \"name=value;...\""
    "--outliers"
        help = "Outlier removal method for training data (none|fence)"
        default = "fence"
    "--mcmc", "-m"
        help = "MCMC samples for hyperparameters"

    # mcmc:
    "--anneal", "-a"
        help = "Annealed samples"
        arg_type = Int
        default = 100
    "--burnin", "-b"
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
    "--silent"
        help = "Don't output progress to stdout (automatic if -o is stdout)"
        action = :store_true
    # predict, sample:
    "--tdata", "-t"
        help = "Data files at which to make predictions"
    "--at"
        help = "Predict at these points (alternative to tdata); format: ;-separated variable=start:step:end OR variable=value"
    # plotfit:
    "--component", "--comp"
        help = "Plot fit plots for components with the given variables"
    # select:
    "--mcmc2", "-m2"
        help = "MCMC samples for hyperparameters of the second model (for model comparison)"

    "--output", "-o"
        help = "Output filename, accepts \"stdout\". For tabular output, supports format flags"
        default = "stdout"

    "mcmc"
        help = "Run MCMC to optimize hyperparameters"
        action = :command
    "predict"
        help = "Calculate the posterior of a GP at tdata given data"
        action = :command
    "sample"
        help = "Sample the posterior of a GP"
        action = :command
    "plotfit"
        help = "Diagnostic plots showing the posteriors of different components of the GP"
        action = :command
    "select"
        help = "Output model selection parameters; requires --mcmc and --mcmc2"
        action = :command
end



function parse_formula(f)
    # TODO 2
end

function read_data(data)
    df = DataFrame()
    key = []
    for file in split(data, ";")
        # (#?[,~]?ID?:)?filename
        m = match(r"^((?<tr>#?)(?<ct>[,~]?)(?<id>[A-Za-z0-9]*):)?(?<file>[^:]*)$", file)

        # Break into parts
        f = m[:file]
        tr = isa(m[:tr], Nothing) ? false : m[:tr] == "#"
        delim = isa(m[:ct], Nothing) ? '\t' : (m[:ct] == "," ? ',' : '\t')
        id = isa(m[:id], Nothing) ? "" : m[:id]

        # Read the data
        # TODO: Errors
        tbl = CSV.read(file, transpose=tr, delim=delim)

        if isempty(df)
            if id == ""
                key = tbl[1]
            else
                # TODO: Errors
                key = tbl[Symbol(id)]
            end
            df = tbl
        else
            if id == ""
                n = zeros(size(tbl)[2])
                for i in 1:size(tbl)[2]
                    n[i] = sum([isa(findfirst(isequal(entry), key), Nothing) for entry in tbl[i]])
                end
                key2 = tbl[findfirst(n .== maximum(n))]
            else
                # TODO: Errors
                key2 = tbl[Symbol(id)]
            end
            ix = [findfirst(isequal(k), key2) for k in key]

            # TODO: Check no 100% match

            tbl2 = tbl[ix,:]

            # TODO: Merge
        end
    end
end

function read_tdata(args)
    # TODO 1: read --tdata or --at
end

function filter_outliers(data, gp, method)
    if method == "none"
        # TODO 3
    elseif method == "fence"
        # TODO 3
    else
        error(@sprintf("Unknown outlier filtering method: %s", method))
    end
end


function read_mcmc(filename)
    # TODO
end

function write_mcmc(chain, filename, append)
    # TODO
end

function read_gp_data(args)
    # Read in input data
    if haskey(args, "data")
        data = read_data(args["data"])
    else
        data = DataFrame()
    end

    # Parse the formula
    gp, data2xyz = parse_formula(args["formula"])

    # Filter outliers
    data = filter_outliers(data, gp, args["outliers"])

    # Transform to GP training vectors x, y, z
    x, y, z = data2xyz(data)
end



function cmd_mcmc(args)
    # Load training data and GP
    gp, x, y, z = read_gp_data(args)

    # Extend a chain?
    if haskey(args, "mcmc")
        chain1 = read_mcmc(args["mcmc"])
        start_θ = chain1[end, gp.chain_names()]
        burnin = 0
        anneal = 0
    else
        start_θ = []
        burnin = args["burnin"]
        anneal = args["anneal"]
    end

    # Run the chain
    chain = gp_mcmc(start_θ, anneal)

    # Thinning
    chain = chain[burnin:thin:end]

    # Output chains
    write_mcmc(chain, args["output"], !isempty(start_θ))
end

function cmd_predict(args)
    # Load training data and GP
    gp, x, y, z = read_gp_data(args)

    # Use fit hyperparameters?
    if haskey(args, "mcmc")
        mcmc = read_mcmc(args["mcmc"])
    else
        mcmc = Chains()
    end

    # Where to sample?
    tdata = read_tdata(args)

    # TODO: Get the predicted latents

    # TODO: Carry the predicted latents through the likelihood

    # TODO: Output predictions
end

function cmd_sample(args)
    # Load training data and GP
    gp, x, y, z = read_gp_data(args)

    # Use fit hyperparameters?
    if haskey(args, "mcmc")
        mcmc = read_mcmc(args["mcmc"])
    else
        mcmc = Chains()
    end

    # Where to sample?
    tdata = read_tdata(args)

    # TODO: Get the predicted latents

    # TODO: Sample the latents

    # TODO: Sample the observation given the latents

    # TODO: Output the sampled values
end

function cmd_plotfit(args)
    # Load training data and GP
    gp, x, y, z = read_gp_data(args)

    # Use fit hyperparameters?
    if haskey(args, "mcmc")
        mcmc = read_mcmc(args["mcmc"])
    else
        mcmc = Chains()
    end

    # For each GP component, output fit plots
    if haskey(args, "component")
        # TODO
    end

    for #each component
        if #this is a component to output
            # TODO make the fitplot for the component
        end
    end
end

function cmd_select(args)
    if !haskey(args, "mcmc")
        error("Expected --mcmc")
    end
    if !haskey(args, "mcmc2")
        error("Expected --mcmc2")
    else

    # Load chains
    mcmc = read_mcmc(args["mcmc"])
    mcmc2 = read_mcmc(args["mcmc2"])

    # Calculate measures
    # TODO

    # Output measures
    # TODO
end
