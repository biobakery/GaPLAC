using Base: Symbol, Real
abstract type GPCompnent end

"""
    SExp(symbol, lengthscale=1)
    SExp(symbol; l = 1)

Squared Exponential Kernel with lengthscale `l`
"""
struct SExp <: GPCompnent
    symbol::Symbol
    lengthscale::Real

    # inner constructor
    SExp(x; l=1) = new(Symbol(x), l)
end

struct OU <: GPCompnent
    symbol::Symbol
    lengthscale::Real

    # inner constructor
    OU(x; l=1) = new(Symbol(x), l)
end

struct Cat <: GPCompnent
    symbol::Symbol
end

Base.+(c1::GPCompnent, c2::GPCompnent) = (:add, c1, c2)
Base.*(c1::GPCompnent, c2::GPCompnent) = (:multiply, c1, c2)