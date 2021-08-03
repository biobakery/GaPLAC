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

function _df_output(df, args)
    if !isnothing(args["output"])
        out = args["output"]
        delim = endswith(out, "csv") ? ',' :
              endswith(out, "tsv") ? '\t' : error("--output arg must be '.tsv' or '.csv'")
        CSV.write(expanduser(args["output"]), df; delim)
    else
        @show df
    end
end