module MCMC

using GaPLAC
using Turing

function run(args)
    @info "running 'mcmc'" 
    gpspec = GaPLAC.gp_formula(args["formula"])
    
    @debug "GP formula" GaPLAC.formula(gpspec)
        
    df = CSV.read(args["data"], DataFrame)
    
    eq, vars = _apply_vars(GaPLAC.formula(gpspec))
    kernels = _walk_kernel(eq)
    length(vars) == length(kernels) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
    inferable = Symbol.(args["infer"])

    y = df[!, GaPLAC.response(gpspec)]
    x = Matrix(df[!, vars])

    @debug "GP equation" eq 
    @debug "Model variables" vars

    @model function inference_engine(Y, X, eq, inferable)
        ℓ ~ TruncatedNormal(0, 5, 0, 30)
        gp, = _apply_vars(eq; hyperparams=Dict(v=> ℓ for v in inferable))
        
        fx ~ AbstractGPs.FiniteGP(GP(gp), RowVecs(X), 0.1)
        Y .~ Normal.(fx, 1)
    end

    m = inference_engine(y, x, GaPLAC.formula(gpspec), inferable)
    
    chain = sample(m, NUTS(0.65), args["samples"], progress=true)
    _df_output(chain, args)

    return chain
end

function _get_hyperparams(infer, kernels, vars)
    hyperparams = Dict{Symbol,Number}()    
end

end