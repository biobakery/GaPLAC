function sample_plot(spec, df)
    gp, vars = GaPLAC.make_gp(spec)
    
    x = df[!, first(vars)]
    y = df[!, GaPLAC.response(spec)]
    fx = AbstractGPs.FiniteGP(gp, x, 0.1)
    
    pgp = posterior(fx, y)
    
    xmin, xmax = extrema(x) .+ (-1, 1)
    xtest = range(xmin, xmax, length = 100)
    ym, yvar = mean_and_var(pgp, xtest) 
    
    fig, ax, l = scatter(x, y, color=:purple, label="samples")
    # srt = sortperm(x)
    # lines!(x[srt],y[srt], color=:purple)
    
    lines!(xtest, ym, color=:dodgerblue, label="mean posterior")
    band!(xtest, ym .- yvar, ym .+ yvar, color=(:dodgerblue, 0.3), label="var posterior")
    
    ax.xlabel = string(first(vars))
    ax.ylabel = string(GaPLAC.response(spec))
    ax.title = "Sample from posterior, x from $(round(xmin, digits=2)) to $(round(xmax, digits=2))"
    axislegend()
    return fig, ax, l
end