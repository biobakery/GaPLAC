module Sample

using GaPLAC
using GaPLAC.AbstractGPs
using GaPLAC.CairoMakie

function run(args)
    @info "running 'sample'" 
    @info args
    gpspec = GaPLAC.gp_spec(args["spec"])
    @debug "GP spec" GaPLAC.formula(gpspec)
        
    gp, vars = GaPLAC.make_gp(gpspec)
    
    @debug "Model variables" vars
    @debug "GP" gp

    atdict = GaPLAC.getatrange(args["at"])
 
    @debug "Inferred ranges" atdict
    @debug "Number of combinations" (*(length.(values(atdict))...))
    
    df = GaPLAC._make_test_df((atdict[v] for v in vars)...; vars)
    X = RowVecs(Matrix(df))
    df[!, GaPLAC.response(gpspec)] = rand(gp(X, 0.1))
    
    GaPLAC._df_output(df, args)

    if !isnothing(args["plot"])
        if length(vars) > 1
            @warn """
                Auto-generated sample plots not available for multivariable models,
                use `--output` to make a table and plot manually
                """
        else
            @info "Plotting output"
            
            fig, ax, l = GaPLAC.sample_plot(gpspec, df)
            save(args["plot"], fig)
        end
    end

    return df, gp
end

end # module