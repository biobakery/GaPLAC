# Covert from cli into KernelFunctions kernels

_default_range(::SqExponentialKernel) = 1:10
_default_range(::LinearKernel)        = -3:0.1:3
_default_range(::CategoricalKernel)   = [1,2,3]

_hyperparam(::SqExponentialKernel) = :lengthscale
_hyperparam(::LinearKernel) = :intercept
_hyperparam(::CategoricalKernel) = :nothing

_convert2kernel(k::SExp) = with_lengthscale(SqExponentialKernel(), k.lengthscale)
_convert2kernel(k::Linear) = LinearKernel(c=k.intercept)
_convert2kernel(::Cat) = CategoricalKernel()

_convert2kernel(gpc::GPCompnent, pos) = _convert2kernel(gpc) ∘ SelectTransform(pos)
_convert2kernel(gpc::Kernel, pos) = _convert2kernel(gpc) ∘ SelectTransform(pos)

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
        return ⊗
    else
        throw(ArgumentError("Operation $op not yet supported"))
    end
end

_convert2eq(op::Symbol, c1::GPCompnent, c2::GPCompnent) = _convert2eq(op, _convert2kernel(c1), _convert2kernel(c2))
_convert2eq(op::Symbol, c1::GPCompnent, c2::Tuple) = _convert2eq(op, _convert2kernel(c1), _convert2eq(c2...))
_convert2eq(op::Symbol, c1::Tuple, c2::GPCompnent) = _convert2eq(op, _convert2eq(c1...), _convert2kernel(c2))

_convert2eq(op::Symbol, c1::Tuple, c2::Tuple) = _convertop(op)(_convert2eq(c1...), _convert2eq(c2...))
_convert2eq(op::Symbol, c1::Kernel, c2::Tuple) = _convertop(op)(c1, _convert2eq(c2...))
_convert2eq(op::Symbol, c1::Tuple, c2::Kernel) = _convertop(op)(_covert2g(c1...), c2)

_convert2eq(op::Symbol, c1::Kernel, c2::Kernel) = _convertop(op)(c1, c2)
_convert2eq(c::GPCompnent) = _convert2kernel(c)

function _apply_vars(formula::Tuple; hyperparams=nothing)
    vars = _formula_pull_varnames(formula)
    ks = _convert2eq(formula...)
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

_apply_vars(formula::GPCompnent) = (_convert2kernel(formula), [_formula_pull_varnames(formula)])