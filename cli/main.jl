
using CSV
using Printf
using Distributions
using Statistics
using ArgParse

function parse_cmdline()
    s = ArgParseSettings()

    @add_arg_table s begin
        #=
        # All commands:
        "formula"
            help = "Gaussian Process formula"
        "--data", "-d"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;~subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            default = "stdin"
        "--bind", "-b"
            help = "Name bindings, format is \"name=value;...\""
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "fence"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
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
        "--atdata", "-t"
            help = "Data files at which to make predictions"
        "--at"
            help = "Predict at these points (alternative to tdata); format: ;-separated variable=start:step:end OR variable=value"
        # fitplot:
        "--component", "--comp"
            help = "Plot fit plots for components with the given variables"
        # select:
        "--mcmc2", "-m2"
            help = "MCMC samples for hyperparameters of the second model (for model comparison)"

        "--output", "-o"
            help = "Output filename, accepts \"stdout\". For tabular output, supports format flags"
            default = "stdout"
            =#

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
    end

    @add_arg_table s["mcmc"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;~subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            default = "stdin"
            required = true
        "--bind", "-b"
            help = "Name bindings, format is \"name=value;...\""
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "fence"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters; if provided, the chain will be extended"
        "--anneal", "-a"
            help = "Annealed samples"
            arg_type = Int
            default = 100
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
        "--silent"
            help = "Don't output progress to stdout (automatic if -o is stdout)"
            action = :store_true
        "--output", "-o"
            help = "Output filename, accepts \"stdout\", supports --data format flags"
            default = "stdout"
    end

    @add_arg_table s["predict"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;~subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            default = "stdin"
            required = true
        "--bind", "-b"
            help = "Name bindings, format is \"name=value;...\""
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "fence"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters"
        "--atdata", "-t"
            help = "Data files providing variable values at which to make predictions"
        "--at"
            help = "Predict at these points (alternative to tdata); format: ;-separated variable=start:step:end OR variable=value"
        "--output", "-o"
            help = "Output filename, accepts \"stdout\", supports --data format flags"
            default = "stdout"
    end

    @add_arg_table s["sample"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;~subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            default = "stdin"
        "--bind", "-b"
            help = "Name bindings, format is \"name=value;...\""
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "fence"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters"
        "--atdata", "-t"
            help = "Data files providing variable values at which to sample the GP"
        "--at"
            help = "Sample at these points (alternative to tdata); format: ;-separated variable=start:step:end OR variable=value"
        "--output", "-o"
            help = "Output filename, accepts \"stdout\", supports --data format flags"
            default = "stdout"
    end

    @add_arg_table s["fitplot"] begin
        "formula"
            help = "Gaussian Process formula"
            required = true
        "--data", "-x"
            help = "Input data files, accepts \"stdin\". ;-separated, use : to provide additional flags which can be combined: \"#:\" transposes the table, \",:\" reads as CSV, \"~:\" reads as TSV (default). Other characters before : give the column/row to join the tables with, e.g. id:data.tsv;~subjectid:subjects.tsv will use the id column of data.tsv and subjectid row of subjects.tsv."
            default = "stdin"
            required = true
        "--bind", "-b"
            help = "Name bindings, format is \"name=value;...\""
        "--rmv_outliers"
            help = "Outlier removal method for training data (none|fence)"
            default = "fence"
        "--outlier_fields"
            help = ";-separated list of additional fields to include in outlier removal"
            default = ""
        "--outlier_ignore"
            help = ";-separated list of fields to ignore for outlier removal"
            default = ""
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters"
        "--component", "--comp"
            help = "Plot fit plots for components with the given variables"
        "--output", "-o"
            help = "Output filename"
            default = "fitplots.pdf"
    end

    @add_arg_table s["select"] begin
        "--mcmc", "-m"
            help = "MCMC samples for hyperparameters"
        "--mcmc2", "-c"
            help = "MCMC samples for hyperparameters of the second model (for model comparison)"
    end

    return parse_args(s)
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
        tbl = CSV.read(file, transpose=tr, delim=delim)
        @info @sprintf("Read %s %s%s", delim=="," ? "CSV" : "TSV", file, tr ? " (transposed)" : "")

        if isempty(df)
            if id == ""
                # Pick the first field for which each value is unique
                goodkeys = [length(unique(tbl[i])) == size(tbl)[1] for i in 1:size(tbl)[2]]
                !any(goodkeys) && error(@sprintf("No suitable id field found for %s", file))
                keyidx = findfirst(goodkeys)
                @info @sprintf("Using %s for sample ids in %s", names(tbl)[keyidx], file)
                key = tbl[keyidx]
            elseif !(Symbol(id) in names(tbl))
                error(@sprintf("%s is not a field in %s", id, file))
            else
                key = tbl[Symbol(id)]
            end
            df = tbl
        else
            if id == ""
                # Pick the field which matches the sample IDs the best
                n = zeros(size(tbl)[2])
                for i in 1:size(tbl)[2]
                    n[i] = sum([isa(findfirst(isequal(entry), key), Nothing) for entry in tbl[i]])
                end
                bestid = indmax(n)
                id = names(tbl)[bestid]
                key2 = tbl[bestid]
            elseif !(Symbol(id) in names(tbl))
                error(@sprintf("%s is not a field in %s", id, file))
            else
                key2 = tbl[Symbol(id)]
            end

            # Match samples
            ix = [findfirst(isequal(k), key2) for k in key]
            if any(isa.(ix, Nothing))
                error(@sprintf("Not all samples mapped to data in %s", file))
            end
            tbl2 = tbl[ix,:]

            # Merge fields
            for i in 1:size(tbl2)[2]
                if names(tbl)[i] != Symbol(id)
                    df[names(tbl)[i]] = tbl2[i]
                end
            end
        end
    end

    return df
end

function read_atdata(args)
    if haskey(args, "atdata") && haskey(args, "at")
        error("--atdata and --at are mutually exclusive")
    end

    # --atdata has the same format as --data
    haskey(args, "atdata") && return read_data(args[:atdata])

    # --at or --atdata must be specified
    !haskey(args, "at") && error("Expected --atdata or --at")

    # Split into name=value pairs
    df = DataFrame()
    for nameexpr in split(args[:at], ";")
        # name = value
        m = match(r"^(?<name>[^=]+?) *= *(?<expr>.*)$", nameexpr)
        isa(m, Nothing) &&
            error("Format for --at is ;-separated \'<name> = <expression>\'")

        # Evaluate
        name = Symbol(m[:name])
        valuefun = genfun(Meta.parse(m[:expr]), names(df))
        value = valuefun([df[i] for i in 1:size(df)[2]]...)

        # What is it?
        if isa(value, Number)
            # Simple value - just assign
            df[name] = value
        elseif isa(value, Vector)
            # Vector - ensure size is correct and assign
            if isempty(df) || length(value) == size(df)[1]
                df[name] = value
            else
                error("The length of %s (%d) does not match the length of earlier values (%d)", name, length(value), size(df)[1])
            end
        elseif isa(value, UnivariateDistribution)
            # Distribution - draw random values
            df[name] = rand.(repeat(value, inner=size(df)[1]))
        else
            error(@sprintf("Unknown value in --at for %s: %s", name, m[:expr]))
        end
    end
end

function atdata_inputs(args, parsedgp)
    # Get the data frame for the target datapoints
    atdata = read_atdata(args, parsedgp)

    # Build GP inputs for the target data
    x, z = gp_inputs(parsedgp, atdata)

    # Make a dataframe of variables used so we can output that for target samples
    index = atdata[union(parsedgp.xvars, parsedgp.zvars)]

    return df, x, z, index
end

function filter_outliers(data, parsedgp, method)
    if method == "none"
        @info "Outlier filtering disabled"
        return data
    end

    # What fields to use?
    outl_fields = parsedgp.xvars
    if haskey(args, "outlier_fields")
        outl_fields = unique(cat(outl_fields, split(args[:outlier_fields], ";")))
    end
    if haskey(args, "outlier_ignore")
        ignore = split(args[:outlier_ignore], ";")
        outl_fields = [x for x in outl_fields if !(x in ignore)]
    end
    outl_df = data[Symbol.(outl_fields)]

    # Select the outlier method
    if method == "fence"
        @info "Outlier filtering using inner fences"
        f = x -> begin
            xnn = x[!isnan(x)]
            q1, q3 = Statistics.quantile(xnn, 0.25), Statistics.quantile(xnn, 0.75)
            iqr = q3 - q1
            (x .> q3 + 1.5 * iqr) .| (x .< q1 - 1.5 * iqr)
        end
    else
        error(@sprintf("Unknown outlier filtering method: %s", method))
    end

    # Which samples are outliers?
    outlier = falses(size(outl_df)[1])
    for i in 1:size(outl_df)[1]
        if isa(outl_df[i], String)
            @info @sprintf("Skipping %s for outlier removal (categorical)", names(outl_df)[i])
        else
            outlier_i = !isnan(outl_df[i]) .& f(outl_df[i])
            outlier .|= outlier_i
            @info @sprintf("%d outliers identified in %s", sum(outlier_i), names(outl_df)[i])
        end
    end

    # Filter outliers
    @info @sprintf("%d outliers removed total", sum(outlier))
    return data[.!outlier,:]
end


function read_mcmc(filename)
    return read_chains(filename)
end

function write_mcmc(chain, filename, append)
    write_chains(chain, filename, append=append)
end

function read_gp_data(args)
    # Read in input data
    data = DataFrame()
    if haskey(args, "data")
        data = read_data(args["data"])
    end

    # Parse the formula
    !("formula" in keys(args)) && error("--formula expected")
    parsedgp = parse_gp_formula(args["formula"])

    # Filter outliers
    data = filter_outliers(data, parsedgp, args["outliers"])

    # Transform to GP training vectors x, y, z
    x, z, y = gp_inputs(parsedgp, data)

    return parsedgp, x, y, z
end

function write_tabular(df, filename)
    # (#?[,~]?ID?:)?filename
    m = match(r"^((?<tr>#?)(?<ct>[,~]?)(?<id>[A-Za-z0-9]*):)?(?<file>[^:]*)$", file)

    # Break into parts
    f = m[:file]
    tr = isa(m[:tr], Nothing) ? false : m[:tr] == "#"
    delim = isa(m[:ct], Nothing) ? '\t' : (m[:ct] == "," ? ',' : '\t')

    # TODO
    tr && error("Transposed output is currently not supported")

    # Write data
    CSV.write(f, df, delim=delim)
end



function cmd_mcmc(args)
    # Load training data and GP
    parsedgp, x, y, z = read_gp_data(args)

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
        # TODO: make a Chains object with only one record which contains the
        # parsed GP hyperparameters
        error("predict from non-fit parameters is NYI")
    end

    # Where to predict?
    tdata, x2, z2, index = atdata_inputs(args, parsedgp)

    # Perform the prediction
    μ_pred, σ2_pred, Q_pred, f_pred, σ2_pred =
        predict(parsedgp.gp, mcmc, x, y, z, x2;
            quantiles=[0.025, 0.05, 0.1, 0.159, 0.5, 0.841, 0.9, 0.95, 0.975])
    prediction = DataFrame(ymu=μ_pred, ystd=sqrt(σ2_pred),
        yQ025=Q_pred[:,1], yQ050=Q_pred[:,2], yQ100=Q_pred[:,3], yQ159=Q_pred[:,4],
        yQ500=Q_pred[:,5],
        yQ841=Q_pred[:,6], yQ900=Q_pred[:,7], yQ950=Q_pred[:,8], yQ975=Q_pred[:,9],
        fmu=f_pred, fstd=sqrt(σ2_pred))

    # Output the sampled values
    write_tabular(hcat(index, prediction), args["output"])
end

function cmd_sample(args)
    # Load training data and GP
    parsedgp, x, y, z = read_gp_data(args)

    # Use fit hyperparameters?
    if haskey(args, "mcmc")
        mcmc = read_mcmc(args["mcmc"])

        # Pick a random set of hyperparameters from the MCMC chain
        ϕ = unrecord(parsedgp.gp, mcmc, rand(DiscreteUniform(1, size(mcmc.df, 1))))
    else
        # Use the parsed parameters
        ϕ = invlink(parsedgp.gp, vcat(parsedgp.θl, parsedgp.θc))
    end

    # Where to sample?
    tdata, x2, z2, index = atdata_inputs(args, parsedgp)

    # Sample the GP
    y2 = samplegp(gp, ϕ, x, y, z, x2, z2)

    # Output the sampled values
    write_tabular(hcat(index, DataFrame(y=y2)), args["output"])
end

function cmd_fitplot(args)
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

    #for each component to output
        #if this is a component to output
            # TODO make the fitplot for the component
        #end
    #end
end

function cmd_select(args)
    # Load chains
    !haskey(args, "mcmc") && error("Expected --mcmc")
    !haskey(args, "mcmc2") && error("Expected --mcmc2")
    mcmc = read_mcmc(args["mcmc"])
    mcmc2 = read_mcmc(args["mcmc2"])

    # Calculate measures
    l2bf = log2_bayes_factor(mcmc, mcmc2)

    # Output to stdout
    println("l2bf")
    println(@sprintf("%.3f", l2bf))
end



# Main
args = parse_cmdline()
cmd = args["%COMMAND%"]
if cmd == "sample"
    cmd_sample(args)
elseif cmd == "predict"
    cmd_predict(args)
elseif cmd == "mcmc"
    cmd_mcmc(args)
elseif cmd == "select"
    cmd_select(args)
elseif cmd == "fitplot"
    cmd_fitplot(args)
end
