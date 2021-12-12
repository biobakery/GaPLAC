# Covert from cli into KernelFunctions kernels

_default_range(::SqExponentialKernel) = 1:10
_default_range(::ExponentialKernel)   = 1:10
_default_range(::LinearKernel)        = -3:0.1:3
_default_range(::CategoricalKernel)   = [1,2,3]

makekernel(k::SExp) = k.lengthscale == 1 ? SqExponentialKernel() : with_lengthscale(SqExponentialKernel(), k.lengthscale)
makekernel(k::OU) = k.lengthscale == 1 ? ExponentialKernel() : with_lengthscale(ExponentialKernel(), k.lengthscale)
makekernel(k::Linear) = LinearKernel(c=k.intercept)
makekernel(::Cat) = CategoricalKernel()

makekernel(::SExp, l) = l == 1 ? SqExponentialKernel() : with_lengthscale(SqExponentialKernel(), l)
makekernel(::OU, l) = l == 1 ? ExponentialKernel() : with_lengthscale(ExponentialKernel(), l)
makekernel(::Linear, c) = LinearKernel(; c)

_walk_kernel(ks::Union{KernelTensorProduct, KernelSum}) = reduce(vcat, [_walk_kernel(k) for k in ks.kernels])
_walk_kernel(ks::TransformedKernel) = _walk_kernel(ks.kernel)
_walk_kernel(k::KernelFunctions.SimpleKernel) = [k]

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
_convert2eq(gpo::GPOperation; hyperparams=Dict()) = _convert2eq(gpo.op, gpo.lhs, gpo.rhs; hyperparams)
_convert2eq(c::GPCompnent; hyperparams=Dict()) = haskey(hyperparams, varname(c)) ? makekernel(c, hyperparams[varname(c)]) : makekernel(c)

_convert2eq(op::Symbol, c1::Kernel, c2::Kernel; hyperparams=Dict()) = _convertop(op)(c1, c2)

# # I think this started from https://julialang.zulipchat.com/#narrow/stream/234072-probprog/topic/Something.20like.20.22Random.20Effects.20Kernel.22.20in.20a.20FiniteGP.3F/near/247353946
# function _apply_select(formula::Tuple;  hyperparams=Dict())
#     vars = varnames(formula)
#     varmap = Dict(var=>i for (i, var) in enumerate(unique(vars)))
    
# end

# formerly _apply_vars
function kernel(formula::Tuple; hyperparams=Dict())
    vars = varnames(formula)
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

kernel(formula::GPCompnent; hyperparams=Dict()) = (_convert2eq(formula; hyperparams), varnames(formula))
