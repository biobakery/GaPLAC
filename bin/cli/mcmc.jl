module MCMC

using GaPLAC
using GaPLAC.AbstractGPs
using GaPLAC.Turing
using GaPLAC.CSV
using GaPLAC.DataFrames
using GaPLAC.Distributions
using Statistics

function run(args)
    @info "running 'mcmc'" 
    gpspec = GaPLAC.gp_spec(args["formula"])
    
    @debug "GP formula" GaPLAC.formula(gpspec)
        
    df = CSV.read(args["data"], DataFrame)
    
    gp, vars = GaPLAC.make_gp(gpspec)
    inferable = Symbol.(args["infer"])
        
    @debug "Model variables" vars
    @debug "GP" gp

    y = df[!, GaPLAC.response(gpspec)]
    x = Matrix(df[!, vars])

    @debug "GP equation" eq 
    @debug "Model variables" vars

    @model function inference_engine(Y, X, eq, inferable)
        ℓ ~ Uniform(0,20)
        gp, = GaPLAC.kernel(eq; hyperparams=Dict(v=> ℓ for v in inferable))
        
        fx ~ AbstractGPs.FiniteGP(GP(gp), RowVecs(X), 0.1)
        Y .~ Normal.(fx, 1)
    end

    m = inference_engine(y, x, GaPLAC.formula(gpspec), inferable)
    
    chain = sample(m, NUTS(0.65), args["samples"], progress=true)
    GaPLAC._df_output(chain, args)

    return chain
end

function _get_hyperparams(infer, kernels, vars)
    hyperparams = Dict{Symbol,Number}()    
end

end