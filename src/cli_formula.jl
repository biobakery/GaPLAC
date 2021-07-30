function _cli_formula_parse(formula::AbstractString)
    spl1 = findfirst(==(':'), formula)
    spl2 = last(findfirst(==('~'), formula))
    isnothing(spl2) && throw(ArgumentError("Invalid formula specification"))
    barind = nextind(formula, spl2)
    formula[barind] == '|' || throw(ArgumentError("Invalid formula specification"))

    if isnothing(spl1)
        lik = GaPLAC.Gaussian()
    else
        spl1 < spl2 || throw(ArgumentError("Invalid formula specification"))
        lik = strip(formula[nextind(formula, spl1):prevind(formula, spl2)])
        lik = isempty(lik) ? GaPLAC.Gaussian() : GaPLAC.eval(Meta.parse(lik))
    end

    response = strip(formula[1:prevind(formula, spl1)])
    
    gp = strip(formula[nextind(formula, barind):end])
    gp = GaPLAC.eval(Meta.parse(gp))

    return Symbol(response), lik, gp
end


function _make_test_df(args...; vars)
    DataFrame(permutedims(reshape(collect(Iterators.flatten(reduce(vcat, Iterators.product(args...)))),
        length(args), *(length.(args)...))), vars)
end

function _cli_run_sample(args)
    @info "running 'sample'" 
    (response, lik, gp_form) = _cli_formula_parse(args["formula"])
    @debug "GP formula" gp_form
        
    eq, vars = _apply_vars(gp_form)
    kernels = _walk_kernel(eq)
    length(vars) == length(kernels) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
    
    @debug "GP equation" eq 
    @debug "Model variables" vars
    gp = GP(eq)

    @debug "GP" gp

    atdict = Dict{Symbol, Any}()
    ats = Meta.parse(args["at"])
    @debug "Passed ranges" ats

    exprs = ats.head == :toplevel ? ats.args : [ats]
    
    for expr in exprs
        expr.head == Symbol("=") || error("Only assignments allowed in `--at` argument")
        var = expr.args[1]
        val = Main.eval(expr.args[2])
        atdict[var] = val
    end
    
    for (i, var) in enumerate(vars)
        if !in(var, keys(atdict))
            atdict[var] = _default_range(kernels[i])
        end
    end
 
    @debug "Inferred ranges" atdict
    @debug "Number of combinations" (*(length.(values(atdict))...))
    
    df = _make_test_df((atdict[v] for v in vars)...; vars)
    X = RowVecs(Matrix(df))
    df[!, response] = rand(gp(X, 0.1))
    
    if !isnothing(args["output"])
        out = args["output"]
        delim = endswith(out, "csv") ? ',' :
              endswith(out, "tsv") ? '\t' : error("--output arg must be '.tsv' or '.csv'")
        CSV.write(expanduser(args["output"]), df; delim)
    else
        @show df
    end

    if !isnothing(args["plot"])
        if length(vars) > 1
            @warn """
                Auto-generated sample plots not available for multivariable models,
                use `--output` to make a table and plot manually
                """
        else
            fig, ax, l = CairoMakie.lines(df[!, first(vars)], df[!, response])
            save(args["plot"], fig)
        end
    end
end

function _cli_run_mcmc(args)
    @info "running 'sample'" 
    (response, lik, gp_form) = _cli_formula_parse(args["formula"])
    @debug "GP formula" gp_form
        
    eq, vars = _apply_vars(gp_form)
    kernels = _walk_kernel(eq)
    length(vars) == length(kernels) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
    
    @debug "GP equation" eq 
    @debug "Model variables" vars
    gp = GP(eq)

end