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
    gpspec = GaPLAC.gp_formula(args["formula"])
    
    @debug "GP formula" GaPLAC.formula(gpspec)
        
    df = CSV.read(args["data"], DataFrame)
    
    eq, vars = GaPLAC._apply_vars(GaPLAC.formula(gpspec))
    kernels = GaPLAC._walk_kernel(eq)
    length(vars) == length(kernels) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
    inferable = Symbol.(args["infer"])

    y = df[!, GaPLAC.response(gpspec)]
    x = Matrix(df[!, vars])

    @debug "GP equation" eq 
    @debug "Model variables" vars

    @model function inference_engine(Y, X, eq, inferable)
        ℓ ~ Uniform(0,20)
        gp, = GaPLAC._apply_vars(eq; hyperparams=Dict(v=> ℓ for v in inferable))
        
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