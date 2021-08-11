using Base: Symbol, Real
abstract type GPCompnent end

struct CategoricalKernel <: KernelFunctions.SimpleKernel end
KernelFunctions.kappa(::CategoricalKernel, d2::Real) = d2 > 0 ? 0.0 : 1.0
KernelFunctions.metric(::CategoricalKernel) = Euclidean()

"""
    SExp(symbol, lengthscale=1)
    SExp(symbol; l = 1)

Squared Exponential Kernel with lengthscale `l`
"""
struct SExp <: GPCompnent
    varname::Symbol
    lengthscale::Real

    # inner constructor
    SExp(x; l=1) = new(Symbol(x), l)
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

Base.:+(c1::GPCompnent, c2::GPCompnent) = (:add, c1, c2)
Base.:+(c1::Tuple, c2::GPCompnent) = (:add, c1, c2)
Base.:+(c1::GPCompnent, c2::Tuple) = (:add, c1, c2)

Base.:*(c1::GPCompnent, c2::GPCompnent) = (:multiply, c1, c2)
Base.:*(c1::GPCompnent, c2::Tuple) = (:multiply, c1, c2)
Base.:*(c1::Tuple, c2::GPCompnent) = (:multiply, c1, c2)
