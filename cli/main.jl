using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Plots
using Logging
using CSV
using Printf
using Distributions
using Statistics
using ArgParse

include("../package/src/gp.jl")
include("../package/src/formula.jl")
include("../package/src/chains.jl")
include("../package/src/mcmc.jl")
include("../package/src/mcmcgp.jl")
include("../package/src/laplacegp.jl")
include("../package/src/directgp.jl")
include("../package/src/bf.jl")

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


function read_data(data)
    df = DataFrame()
    key = []
    for file in split(data, ";")
        # (#?[,~]?ID?:)?filename
        m = match(r"^((?<tr>#?)(?<ct>[,~]?)(?<id>[A-Za-z0-9]*):)?(?<file>[^:]*)$", file)

        # Break into parts
        filename = m[:file]
        tr = isa(m[:tr], Nothing) ? false : m[:tr] == "#"
        delim = if isa(m[:ct], Nothing)
            endswith(filename, "csv") ? ',' : '\t'
        else
            m[:ct] == "," ? ',' : '\t'
        end
        id = isa(m[:id], Nothing) ? "" : m[:id]

        # Read the data
        file = filename == "stdin" ? stdin : string(filename)
        tbl = CSV.read(file, transpose=tr, delim=delim)
        @info @sprintf("Read %s %s%s (%d samples with %d features)",
            delim=="," ? "CSV" : "TSV", file, tr ? " (transposed)" : "",
            size(tbl,1), size(tbl,2))

        if isempty(df)
            if id == ""
                # Pick the first field for which each value is unique
                goodkeys = [length(unique(tbl[i])) == size(tbl,1) for i in 1:size(tbl,2)]
                !any(goodkeys) && error(@sprintf("No suitable id field found for %s", file))
                keyidx = findfirst(goodkeys)
                @info @sprintf("Using %s for sample ids in %s", names(tbl)[keyidx], file)
                key = tbl[:,keyidx]
            elseif !(Symbol(id) in names(tbl))
                error(@sprintf("%s is not a field in %s", id, file))
            else
                key = tbl[:,Symbol(id)]
            end
            df = tbl
        else
            if id == ""
                # Pick the field which matches the sample IDs the best
                n = zeros(size(tbl, 2))
                for i in 1:size(tbl, 2)
                    n[i] = sum([isa(findfirst(isequal(entry), key), Nothing) for entry in tbl[i]])
                end
                bestid = indmax(n)
                id = names(tbl)[bestid]
                key2 = tbl[:,bestid]
            elseif !(Symbol(id) in names(tbl))
                error(@sprintf("%s is not a field in %s", id, file))
            else
                key2 = tbl[:,Symbol(id)]
            end

            # Match samples
            ix = [findfirst(isequal(k), key2) for k in key]
            if any(i-> isa(i, Nothing), ix)
                error(@sprintf("Not all samples mapped to data in %s", file))
            end
            tbl2 = tbl[ix,:]

            # Merge fields
            for i in 1:size(tbl2, 2)
                if names(tbl)[i] != Symbol(id)
                    df[!, names(tbl)[i]] = tbl2[:,i]
                end
            end
        end
    end

    return df
end

function read_atdata(args)
    if !isa(args["atdata"], Nothing) && !isa(args["at"], Nothing)
        error("--atdata and --at are mutually exclusive")
    end

    # --atdata has the same format as --data
    !isa(args["atdata"], Nothing) && return read_data(args[:atdata])

    # --at or --atdata must be specified
    isa(args["at"], Nothing) && error("Expected --atdata or --at")

    # Split into name=value pairs
    df = DataFrame()
    for nameexpr in split(args["at"], ";")
        # name = value
        m = match(r"^(?<name>[^=]+?) *(/(?<group>[^=]+))? *= *(?<expr>.*)$", nameexpr)
        isa(m, Nothing) &&
            error("Format for --at is ;-separated \'<name> = <expression>\'")

        # Evaluate
        name = Symbol(m[:name])
        valuefun = genfun(Meta.parse(m[:expr]), names(df))
        value = Base.invokelatest(valuefun, [df[:,i] for i in 1:size(df)[2]]...)

        # What is it?
        if isa(value, Number)
            # Simple value - just assign
            df[name] = value
        elseif isa(value, Vector) || isa(value, UnitRange) || isa(value, StepRange) || isa(value, StepRangeLen)
            # Vector - ensure size is correct and assign
            value = collect(value)
            if !isa(m[:group], Nothing)
                isempty(df) && error("--at cannot have a grouped variable as the first variable")
                groupfun = genfun(Meta.parse(m[:group]), names(df))
                group = Base.invokelatest(groupfun, [df[i] for i in 1:size(df)[2]]...)
                unq_group = unique(group)
                if length(unq_group) == length(group)
                    # Unique groups - expand
                    df = repeat(df, inner=length(value))
                    df[!,name] = repeat(value, outer=length(group))
                else
                    # Non-unique groups - fill in
                    n = sum(group .== unq_group[1])
                    values = repeat(value[1], inner=size(df)[1])
                    n != length(value) && error("Values must have the same length as the groups")
                    for ug in unq_group
                        mask = group .== ug
                        n != sum(mask) && error(@sprintf("Grouping variable %s has unequal group lengths", m[:group]))
                        values[mask] = value
                    end
                    df[name] = values
                end
            elseif isempty(df) || length(value) == size(df, 1)
                df[!, name] = value
            else
                error(@sprintf("The length of %s (%d) does not match the length of earlier values (%d)", name, length(value), size(df)[1]))
            end
        elseif isa(value, UnivariateDistribution)
            # Distribution - draw random values
            df[name] = rand.(repeat([value], inner=size(df)[1]))
        else
            error(@sprintf("Unknown value in --at for %s: %s", name, m[:expr]))
        end
    end

    return df
end

function filter_outliers(data, parsedgp, args)
    if size(data)[1] == 0
        # If there's no data, don't log anything related to outlier filtering
        return data
    end

    # How to filter outliers
    method = args["rmv_outliers"]
    if method == "none"
        @info "Outlier filtering disabled"
        return data
    end

    # What fields to use?
    outl_fields = string.([union(parsedgp.xfun_params, parsedgp.yfun_params)...])
    if !isa(args["outlier_fields"], Nothing)
        outl_fields = unique(vcat(outl_fields, split(args["outlier_fields"], ";")))
    end
    if !isa(args["outlier_ignore"], Nothing)
        ignore = split(args["outlier_ignore"], ";")
        outl_fields = [x for x in outl_fields if !(x in ignore)]
    end
    outl_df = data[Symbol.(outl_fields)]

    # Select the outlier method
    if method == "fence"
        @info "Outlier filtering using inner fences"
        f = x -> begin
            xnn = x[.!isnan.(x)]
            q1, q3 = Statistics.quantile(xnn, 0.25), Statistics.quantile(xnn, 0.75)
            iqr = q3 - q1
            (x .> q3 + 1.5 * iqr) .| (x .< q1 - 1.5 * iqr)
        end
    else
        error(@sprintf("Unknown outlier filtering method: %s", method))
    end

    # Which samples are outliers?
    outlier = falses(size(outl_df)[1])
    for i in 1:size(outl_df)[2]
        if isa(outl_df[i], String)
            @info @sprintf("Skipping %s for outlier removal (categorical)", names(outl_df)[i])
        else
            outlier_i = .!isnan.(outl_df[i]) .& f(outl_df[i])
            outlier .|= outlier_i
            @info @sprintf("%d outliers identified in %s", sum(outlier_i), names(outl_df)[i])
        end
    end

    # Filter outliers
    @info @sprintf("%d outliers removed total", sum(outlier))
    return data[.!outlier,:]
end


function read_mcmc(filename)
    file = filename == "stdin" ? stdin : filename
    chain = read_chains(file)

    @info @sprintf("Read %d MCMC samples from %s", size(chain.df,1), filename)
    return chain
end

function write_mcmc(chain, filename, append)
    file = filename == "stdout" ? stdout : filename
    write_chains(chain, file; append=append)

    @info @sprintf("%s %d MCMC samples to %s", append ? "Appended" : "Wrote", size(chain.df,1), filename)
end

function read_gp_data(args, need_at=false)
    # Read in input data
    data = DataFrame()
    if !isa(args["data"], Nothing)
        data = read_data(args["data"])
    else
        @info "No training points specified"
    end
    atdata = DataFrame()
    if need_at
        atdata = read_atdata(args)
    end

    # Parse the formula
    isa(args["formula"], Nothing) && error("--formula expected")
    @info @sprintf("Formula: %s", args["formula"])
    varnames = []
    iscat = []
    if !isempty(data)
        varnames = string.(names(data))
        iscat = [!(typeof(data[i][1]) <: Number) for i in 1:size(data,2)]
    elseif !isempty(atdata)
        varnames = string.(names(atdata))
        iscat = [!(typeof(atdata[i][1]) <: Number) for i in 1:size(atdata,2)]
    end
    parsedgp = parse_gp_formula(args["formula"], varnames, iscat)

    # Filter outliers
    data = filter_outliers(data, parsedgp, args)

    return parsedgp, data, atdata
end

function gen_gp_inputs(parsedgp, data)
    # Transform to GP training vectors x, y, z
    x, z, y = gp_inputs(parsedgp, data)
    return x, y, z
end

function gen_gp_inputs(parsedgp, data, atdata)
    # Transform to GP training vectors x, y, z
    x, z, y = gp_inputs(parsedgp, data)

    # Get target data matrices
    x2, z2 = gp_inputs(parsedgp, atdata)
    vars = union(parsedgp.xfun_params, parsedgp.zfun_params)
    index = atdata[:,[(name in vars) for name in names(atdata)]] # returns a copy

    return x, y, z, x2, z2, index
end

function write_tabular(df, filespec)
    # (#?[,~]?ID?:)?filename
    m = match(r"^((?<tr>#?)(?<ct>[,~]?)(?<id>.*):)?(?<file>[^:]*)$", filespec)

    # Break into parts
    filename = m[:file]
    tr = isa(m[:tr], Nothing) ? false : m[:tr] == "#"
    delim = isa(m[:ct], Nothing) ? '\t' : (m[:ct] == "," ? ',' : '\t')

    # TODO
    tr && error("Transposed output is currently not supported")

    # Write data
    if filename == "stdout"
        # CSV.write calls seek on the file, so for stdout which doesn't
        # support seed, we need a different call
        CSV.write(stdout, df, delim=delim, append=true, writeheader=true)
    else
        CSV.write(string(filename), df, delim=delim)
    end
end

function make_sampleplot(args, data, y, parsedgp)
    # Get the time and group vectors based on --plotx
    time, group, xlab = [], [], "x"
    if isa(args["plotx"], Nothing)
        nunq = [length(unique(data[i])) for i in 1:size(data)[2]]
        ni = findfirst(nunq.==maximum(nunq))
        time = data[ni]
        group = zeros(length(time))
        xlab = names(data)[ni]
    else
        m = match(r"^(?<time>.*?) *(: *(?<group>.*))?$", args["plotx"])
        !(string(m[:time]) in string.(names(data))) && error("X dimension of plot is not a variable given to --at or --atdata")
        time = data[!,Symbol(m[:time])]
        group = if isa(m[:group], Nothing)
            zeros(length(time))
        else
            !(string(m[:group]) in string.(names(data))) && error("Grouping variable for plot is not a variable given to --at or --atdata")
            data[!, Symbol(m[:group])]
        end
        xlab = string(m[:time])
    end

    # Make the plot!
    plt = plot(size=(600,400), legend=false, xlabel=xlab, ylabel=parsedgp.yname)
    unq_group = unique(group)
    for ug in unq_group
        # Separate line for each group
        mask = group .== ug
        tg = time[mask]
        yg = y[mask]
        I = sortperm(tg)
        plot!(plt, tg[I], yg[I])
    end

    # Dump to file
    savefig(plt, args["plot"])
end



function cmd_mcmc(args)
    # Load training data and GP
    parsedgp, data = read_gp_data(args)

    # Refresh the world age before continuing command processing
    Base.invokelatest(cmd_mcmc_cont, args, parsedgp, data)
end

function cmd_mcmc_cont(args, parsedgp, data)
    # Get GP inputs
    x, y, z = gen_gp_inputs(parsedgp, data)

    # Extend a chain?
    burnin = 0
    if !isa(args["mcmc"], Nothing)
        chain1 = read_mcmc(args["mcmc"])
        start_θ = chain1[end, gp.chain_names()]
    else
        start_θ = vcat(parsedgp.θc, parsedgp.θl)
        burnin = args["burnin"]
    end

    # Get the properties of the chain
    thinning = args["thin"]
    samples = args["samples"]

    # Run the chain
    chain = mcmcgp(parsedgp.gp, x, y, z, start_θ, burnin + 1 + thinning * (samples-1))

    # Thinning
    chain = thin(chain, burnin, thinning)

    # Output chains
    write_mcmc(chain, args["output"], !isa(args["mcmc"], Nothing))
end

function cmd_predict(args)
    # Load training data and GP
    parsedgp, data, atdata = read_gp_data(args, true)

    # Refresh the world age before continuing command processing
    Base.invokelatest(cmd_predict_cont, args, parsedgp, data, atdata)
end

function cmd_predict_cont(args, parsedgp, data, atdata)
    # Get GP inputs
    x, y, z, x2, z2, index = gen_gp_inputs(parsedgp, data, atdata)

    # Use fit hyperparameters?
    if !isa(args["mcmc"], Nothing)
        mcmc = read_mcmc(args["mcmc"])
    else
        θ = vcat(parsedgp.θc, parsedgp.θl)
        mcmc = Chains()
        record!(parsedgp.gp, mcmc, invlink(parsedgp.gp, θ))
    end

    # Perform the prediction
    μ_pred, σ2_pred, Q_pred, f_pred, σ2_pred =
        predict(parsedgp.gp, mcmc, x, y, z, x2, z2;
            quantiles=[0.025, 0.05, 0.1, 0.159, 0.5, 0.841, 0.9, 0.95, 0.975])
    prediction = DataFrame(ymu=μ_pred, ystd=sqrt.(σ2_pred),
        yQ025=Q_pred[:,1], yQ050=Q_pred[:,2], yQ100=Q_pred[:,3], yQ159=Q_pred[:,4],
        yQ500=Q_pred[:,5],
        yQ841=Q_pred[:,6], yQ900=Q_pred[:,7], yQ950=Q_pred[:,8], yQ975=Q_pred[:,9],
        fmu=f_pred, fstd=sqrt.(σ2_pred))

    # Output the sampled values
    write_tabular(hcat(index, prediction), args["output"])
end

function cmd_sample(args)
    # Load training data and GP
    parsedgp, data, atdata = read_gp_data(args, true)

    # Refresh the world age before continuing command processing
    Base.invokelatest(cmd_sample_cont, args, parsedgp, data, atdata)
end

function cmd_sample_cont(args, parsedgp, data, atdata)
    # Get GP inputs
    x, y, z, x2, z2, index = gen_gp_inputs(parsedgp, data, atdata)

    # Use fit hyperparameters?
    if !isa(args["mcmc"], Nothing)
        mcmc = read_mcmc(args["mcmc"])

        # Pick a random set of hyperparameters from the MCMC chain
        ϕ = unrecord(parsedgp.gp, mcmc, rand(DiscreteUniform(1, length(mcmc))))
    else
        # Use the parsed parameters
        ϕ = invlink(parsedgp.gp, vcat(parsedgp.θl, parsedgp.θc))
    end

    # Sample the GP
    y2, f2 = samplegp(parsedgp.gp, ϕ, x, y, z, x2, z2)

    # Output a plot
    if !isa(args["plot"], Nothing)
        make_sampleplot(args, atdata, y2, parsedgp)
    end

    # Output the sampled values
    write_tabular(hcat(index, DataFrame(y=y2, f=f2)), args["output"])
end

function cmd_fitplot(args)
    # Load training data and GP
    parsedgp, data = read_gp_data(args)

    # Refresh the world age before continuing command processing
    Base.invokelatest(cmd_fitplot_cont, args, parsedgp, data)
end

function cmd_fitplot_cont(args, parsedgp, data)
    # Get GP inputs
    x, y, z = gen_gp_inputs(parsedgp, data)

    # Use fit hyperparameters?
    if !isa(args["mcmc"], Nothing)
        mcmc = read_mcmc(args["mcmc"])
    else
        mcmc = Chains()
    end

    # For each GP component, output fit plots
    if !isa(args["component"], Nothing)
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
    isa(args["mcmc"], Nothing) && error("Expected --mcmc")
    isa(args["mcmc2"], Nothing) && error("Expected --mcmc2")
    mcmc = read_mcmc(args["mcmc"])
    mcmc2 = read_mcmc(args["mcmc2"])

    # Calculate measures
    l2bf = log2_bayes_factor(mcmc, mcmc2)

    # Output to stdout
    println(@sprintf("%.3f", l2bf))
end

function processcmd(args)
    # Process commands
    cmd = args["%COMMAND%"]
    if cmd == "sample"
        cmd_sample(args[cmd])
    elseif cmd == "predict"
        cmd_predict(args[cmd])
    elseif cmd == "mcmc"
        cmd_mcmc(args[cmd])
    elseif cmd == "select"
        cmd_select(args[cmd])
    elseif cmd == "fitplot"
        cmd_fitplot(args[cmd])
    end
end


# Main
if !Base.isinteractive()
    args = parse_cmdline()

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
