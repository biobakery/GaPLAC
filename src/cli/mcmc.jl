function _cli_run_mcmc(args)
    @info "running 'mcmc'" 
    (response, lik, gp_form) = _cli_formula_parse(args["formula"])
    @debug "GP formula" gp_form
        
    df = CSV.read(args["data"], DataFrame)
    
    eq, vars = _apply_vars(gp_form)
    kernels = _walk_kernel(eq)
    length(vars) == length(kernels) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
    
    y = df[!, response]
    x = Matrix(df[!, vars])

    @debug "GP equation" eq 
    @debug "Model variables" vars


    @model function inference_engine(Y, X, ks)
        
        k = σ * with_lengthscale(first(ks), ℓ)
        
        gp = GP(k)
        fx ~ AbstractGPs.FiniteGP(gp, RowVecs(X), σ)    
        
        Y .~ Normal.(fx, σ)
    end

    m = inference_engine(x, y, eq)
    
    chain = sample(m, NUTS(0.65), 200, progress=true)
    _df_output(chain, args)

    @info chain
end