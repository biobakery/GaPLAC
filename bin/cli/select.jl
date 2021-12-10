module Select

using GaPLAC
using GaPLAC.CSV
using GaPLAC.DataFrames
using GaPLAC.StatsBase
using GaPLAC.AbstractGPs

function run(args)
    @info "running 'select'"
    @info args 
    
    !isempty(args["chains"]) && !isempty(args["formulae"]) && throw(ArgumentError("'select' can only take one of '--formulae' or '--chains', not both"))

    if !isempty(args["chains"])
        lp1df = CSV.read(args["chains"][1], DataFrame)
        lp1 = log2(harmmean([BigFloat(2) ^ x for x in lp1df[!, :lp]]))
        lp2df = CSV.read(args["chains"][2], DataFrame)
        lp2 = log2(harmmean([BigFloat(2) ^ x for x in lp2df[!, :lp]]))
        bayes = log2(BigFloat(2)^lp1 / BigFloat(2)^lp2)
    elseif !isempty(args["formulae"])
        spec1 = GaPLAC.gp_formula(args["formulae"][1])
        spec2 = GaPLAC.gp_formula(args["formulae"][2])
        
        @debug "GP formulae" GaPLAC.formula(spec1) GaPLAC.formula(spec2)
            
        eq1, vars1 = GaPLAC._apply_vars(GaPLAC.formula(spec1))
        kernels1 = GaPLAC._walk_kernel(eq1)
        length(vars1) == length(kernels1) || error("Something went wrong with equation parsing, number of variables should == number of kernels")

        eq2, vars2 = GaPLAC._apply_vars(GaPLAC.formula(spec2))
        kernels2 = GaPLAC._walk_kernel(eq2)
        length(vars2) == length(kernels2) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
        
        gp1 = GP(eq1)
        gp2 = GP(eq2)

        @debug "GPs:" gp1 gp2

        df = CSV.read(args["data"], DataFrame)
        df = disallowmissing(df[completecases(df),:])
        
        y1 = df[!, GaPLAC.response(spec1)]
        x1 = Matrix(df[!, vars1])
        pr1 = AbstractGPs.FiniteGP(gp1, x1, 0.1, obsdim=1)    

        y2 = df[!, GaPLAC.response(spec2)]
        x2 = Matrix(df[!, vars2])
        pr2 = AbstractGPs.FiniteGP(gp2, x2, 0.1, obsdim=1)    
        
        lp1 = logpdf(pr1, y1)
        lp2 = logpdf(pr2, y2)
        pst1 = posterior(pr1, y1)
        pst2 = posterior(pr2, y2)
        
        bayes = log2(BigFloat(2)^lp1 / BigFloat(2)^lp2)
    else
        throw(ArgumentError("'select' command requires either '--chains' or '--formulae' arguments"))
    end

    @info """
        **Log2 Bayes**: $(round(Float64(bayes), digits=3))
        
        - **Log(pdf)** - model 1: $(round(lp1, digits=4))
        - **Log(pdf)** - model 2: $(round(lp2, digits=4))
        
        _Note_ - Positive values indicate more evidence for model 1
        """

end

end