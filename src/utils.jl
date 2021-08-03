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

