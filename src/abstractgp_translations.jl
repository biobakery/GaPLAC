# Covert from cli into KernelFunctions kernels

_default_range(::SqExponentialKernel) = 1:10
_default_range(::LinearKernel)        = -3:0.1:3
_default_range(::CategoricalKernel)   = [1,2,3]


_convert2kernel(k::SExp) = k.lengthscale == 1 ? SqExponentialKernel() : with_lengthscale(SqExponentialKernel(), k.lengthscale)
_convert2kernel(k::OU) = k.lengthscale == 1 ? ExponentialKernel() : with_lengthscale(ExponentialKernel(), k.lengthscale)
_convert2kernel(k::Linear) = LinearKernel(c=k.intercept)
_convert2kernel(::Cat) = CategoricalKernel()

_convert2kernel(::SExp, l) = l == 1 ? SqExponentialKernel() : with_lengthscale(SqExponentialKernel(), l)
_convert2kernel(::OU, l) = l == 1 ? ExponentialKernel() : with_lengthscale(ExponentialKernel(), l)
_convert2kernel(::Linear, c) = LinearKernel(; c)

_walk_kernel(ks::Union{KernelTensorProduct, KernelSum}) = reduce(vcat, [_walk_kernel(k) for k in ks.kernels])
_walk_kernel(ks::TransformedKernel) = _walk_kernel(ks.kernel)
_walk_kernel(k::KernelFunctions.SimpleKernel) = [k]

function _formula_pull_varnames(c::GPCompnent)
    return varname(c)
end

function _formula_pull_varnames(c::Tuple)
    reduce(vcat, _formula_pull_varnames.(c[2:3]))
end

function _convertop(op::Symbol)
    if op == :add
        return +
    elseif op == :multiply
        return *
    else
        throw(ArgumentError("Operation $op not yet supported"))
    end
end

_convert2eq(op::Symbol, c1, c2; hyperparams=Dict()) = _convert2eq(op, _convert2eq(c1; hyperparams), _convert2eq(c2; hyperparams))
_convert2eq(t::Tuple; hyperparams=Dict()) = _convert2eq(t...; hyperparams)
_convert2eq(c::GPCompnent; hyperparams=Dict()) = haskey(hyperparams, varname(c)) ? _convert2kernel(c, hyperparams[varname(c)]) : _convert2kernel(c)

_convert2eq(op::Symbol, c1::Kernel, c2::Kernel; hyperparams=Dict()) = _convertop(op)(c1, c2)

function _apply_select(formula::Tuple;  hyperparams=Dict())
    vars = _formula_pull_varnames(formula)
    varmap = Dict(var=>i for (i, var) in enumerate(unique(vars)))
    
end

function _apply_vars(formula::Tuple; hyperparams=Dict())
    vars = _formula_pull_varnames(formula)
    ks = _convert2eq(formula; hyperparams)
    retkernel = nothing
    counter = 1
    current_k = 1
    while current_k <= length(vars)
        k = ks.kernels[counter]
        if k isa KernelTensorProduct
            kernel = k ∘ SelectTransform(current_k:(current_k + length(k.kernels) - 1))
            retkernel = isnothing(retkernel) ? kernel : retkernel + kernel
            current_k += length(k.kernels) 
            counter += 1
        elseif k isa KernelFunctions.SimpleKernel || k isa KernelFunctions.TransformedKernel
            kernel = k ∘ SelectTransform([current_k])
            retkernel = isnothing(retkernel) ? kernel : retkernel + kernel
            current_k += 1
            counter += 1
        else
            @show typeof(k)
            error()
        end
    end
    return retkernel, vars
end

_apply_vars(formula::GPCompnent; hyperparams=Dict()) = (_convert2eq(formula; hyperparams), [_formula_pull_varnames(formula)])
