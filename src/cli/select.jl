function _cli_run_select(args)
    @info "running 'select'" 
    
    haskey(args, "chains") && haskey(args, "formulae") && throw(ArgumentError("'select' can only take one of '--formulae' or '--chains', not both"))

    if haskey(args, "chains")
        lp1 = CSV.read(args["chains"][1], DataFrame)
        lp2 = CSV.read(args["chains"][2], DataFrame)

        @info log2(harmmean(BigFloat(2) ^ x for x in lp1[!, :lp]) /
                   harmmean(BigFloat(2) ^ x for x in lp2[!, :lp]))
    elseif haskey(args, "formulae")
        (response1, lik1, gp_form1) = _cli_formula_parse(args["formulae"][1])
        (response2, lik2, gp_form2) = _cli_formula_parse(args["formulae"][2])
        
        @debug "GP formulae" gp_form1 gp_form2
            
        eq1, vars1 = _apply_vars(gp_form1)
        kernels1 = _walk_kernel(eq1)
        length(vars1) == length(kernels1) || error("Something went wrong with equation parsing, number of variables should == number of kernels")

        eq2, vars2 = _apply_vars(gp_form2)
        kernels2 = _walk_kernel(eq2)
        length(vars2) == length(kernels2) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
        
        gp1 = GP(eq1)
        gp2 = GP(eq2)

        @debug "GPs:" gp1 gp2

        df = CSV.read(args["data"], DataFrame)

        df1 = unique(df, vars1)
        df1 = disallowmissing(df1[completecases(df1),:])
        y1 = df1[!, response1]
        x1 = Matrix(df1[!, vars1])
        pr1 = AbstractGPs.FiniteGP(gp1, x1, 0.1, obsdim=1)    

        df2 = unique(df, vars2)
        df2 = disallowmissing(df2[completecases(df2),:])
        y2 = df2[!, response2]
        x2 = Matrix(df2[!, vars2])
        pr2 = AbstractGPs.FiniteGP(gp2, x2, 0.1, obsdim=1)    
        
        pst1 = posterior(pr1, y1)
        pst2 = posterior(pr2, y2)
        
        @info "Model 1 / Model 2" log2(BigFloat(2)^logpdf(pr2, y2) / BigFloat(2)^logpdf(pr1, y1))
    else
        throw(ArgumentError("'select' command requires either '--chains' or 'formulae' arguments"))
    end
end