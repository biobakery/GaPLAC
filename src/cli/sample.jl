function _cli_run_sample(args)
    @info "running 'sample'" 
    @info args
    (resp, lik, gp_form) = _cli_formula_parse(args["formula"])
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
        val = GaPLAC.eval(expr.args[2])
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
    df[!, resp] = rand(gp(X, 0.1))
    
    _df_output(df, args)

    if !isnothing(args["plot"])
        if length(vars) > 1
            @warn """
                Auto-generated sample plots not available for multivariable models,
                use `--output` to make a table and plot manually
                """
        else
            @info "Plotting output"
            x = df[!, first(vars)]
            y = df[!, resp]
            fx = AbstractGPs.FiniteGP(gp, x, 0.1)
            pgp = posterior(fx, y)
            xmin, xmax = extrema(x) .+ (-1, 1)
            xtest = range(xmin, xmax, length = 100)
            ym, yvar = mean_and_var(pgp, xtest) 
            
            fig, ax, l = scatter(x, y, color=:purple)
            srt = sortperm(x)
            lines!(x[srt],y[srt], color=:purple)
            ax.xlabel = string(first(vars))
            ax.ylabel = string(resp)
            ax.title = "Sample from posterior, x from $xmin to $xmax"
            lines!(xtest, ym, color=:dodgerblue)
            band!(xtest, ym .- yvar, ym .+ yvar, color=(:dodgerblue, 0.3))

            save(args["plot"], fig)
        end
    end

    return df, gp
end
