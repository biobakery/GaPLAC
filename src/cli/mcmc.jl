function _cli_run_mcmc(args)
    @info "running 'mcmc'" 
    (response, lik, gp_form) = _cli_formula_parse(args["formula"])
    @debug "GP formula" gp_form
        
    df = CSV.read(args["data"], DataFrame)
    
    eq, vars = _apply_vars(gp_form)
    kernels = _walk_kernel(eq)
    length(vars) == length(kernels) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
    inferable = Symbol.(args["infer"])

    y = df[!, response]
    x = Matrix(df[!, vars])

    @debug "GP equation" eq 
    @debug "Model variables" vars


    @model function inference_engine(Y, X, eq, inferable)
        ℓ ~ TruncatedNormal(1, 5, 0, 30)
        gp, = _apply_vars(eq; hyperparams=Dict(v=> ℓ for v in inferable))
        
        fx ~ AbstractGPs.FiniteGP(GP(gp), RowVecs(X), 0.1)
        Y .~ Normal.(fx, 1)
    end

    m = inference_engine(y, x, gp_form, inferable)
    
    chain = sample(m, NUTS(0.65), 200, progress=true)
    _df_output(chain, args)

    @info chain
end

function _get_hyperparams(infer, kernels, vars)
    hyperparams = Dict{Symbol,Number}()    
end