function _make_test_grid(args...)
    return permutedims(reshape(collect(Iterators.flatten(reduce(vcat, Iterators.product(args...)))),
    length(args), *(length.(args)...)))
end

function _make_test_df(args...; vars)
    DataFrame(_make_test_grid(args...), vars)
end

function _make_test_rows(args...)
    RowVecs(_make_test_grid(args...))
end

function getrank(v; flattenzeros=true)
    r = invperm(sortperm(v))
    if flattenzeros
        z = findall(iszero, v)
        r[z] .= 1
    end
    return r
end

function invnormaltransform(v; μ=0, σ=1, c=3/8, flattenzeros=true)
    rank = getrank(v; flattenzeros)
    return [norminvcdf(μ, σ, (x - c) / (length(v) - 2c + 1)) for x in rank]
end

function _df_output(df, args)
    if !isnothing(args["output"])
        @info "Writing to $(args["output"])"
        out = args["output"]
        delim = endswith(out, "csv") ? ',' :
              endswith(out, "tsv") ? '\t' : error("--output arg must be '.tsv' or '.csv'")
        CSV.write(expanduser(args["output"]), df; delim)
    else
        @show df
    end
end
