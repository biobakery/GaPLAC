using Base: Symbol, Real

abstract type GPCompnent end

struct GPOperation <: GPCompnent
    op::Symbol
    lhs::GPCompnent
    rhs::GPCompnent
end

struct CategoricalKernel <: KernelFunctions.SimpleKernel end
KernelFunctions.kappa(::CategoricalKernel, d2::Real) = d2 > 0 ? 0.0 : 1.0
KernelFunctions.metric(::CategoricalKernel) = Euclidean()

"""
    SqExp(symbol, lengthscale=1)
    SqExp(symbol; l = 1)

Squared Exponential Kernel with lengthscale `l`
"""
struct SqExp <: GPCompnent
    varname::Symbol
    lengthscale::Real

    # inner constructor
    SqExp(x; l=1) = new(Symbol(x), l)
end

struct Linear <: GPCompnent
    varname::Symbol
    intercept::Real

    # inner constructor
    Linear(x; c=0) = new(Symbol(x), c)
end

struct OU <: GPCompnent
    varname::Symbol
    lengthscale::Real

    # inner constructor
    OU(x; l=1) = new(Symbol(x), l)
end

struct Cat <: GPCompnent
    varname::Symbol
end

## Accessor Functions

varname(gpc::GPCompnent) = gpc.varname
varnames(gpc::GPCompnent) = [varname(gpc)]
varnames(gpo::GPOperation) = reduce(hcat, varname.(gpo.lhs, gpo.rhs))

Base.:+(c1::GPCompnent, c2::GPCompnent) = GPOperation(:add, c1, c2)
# Base.:+(c1::Tuple, c2::GPCompnent) = (:add, c1, c2)
# Base.:+(c1::GPCompnent, c2::Tuple) = (:add, c1, c2)

Base.:*(c1::GPCompnent, c2::GPCompnent) = GPOperation(:multiply, c1, c2)
# Base.:*(c1::GPCompnent, c2::Tuple) = (:multiply, c1, c2)
# Base.:*(c1::Tuple, c2::GPCompnent) = (:multiply, c1, c2)
