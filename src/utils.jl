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


function invnormaltransform(v; μ=0, σ=1, c=3/8)
    rank = invperm(sortperm(v))
    N = length(v)
    return [norminvcdf(μ, σ, (x - c) / (N - 2c + 1)) for x in rank]
end