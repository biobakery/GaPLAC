
# round(y * filtered_reads) : BetaBinomialR(filtered_reads) ~|
#   Cat(subject) * OU(time; l=Gamma(1)) +
#   Linear(diet)

using DataFrames
using Distributions
using Printf
using SyntaxTree
using ExpressionUtils

struct ParsedGPFormula
    # Function to produce the GP's output
    # Signature: (<yfun_params>) -> Number
    yfun::Function
    # Parameters of the y function
    yfun_params::Array{Symbol}

    # Function to produce the GP's covariance function input matrix
    # Signature: (<xfun_params>) -> Vector{Float64}
    xfun::Function
    # Parameters of the x function
    xfun_params::Array{Symbol}

    # Function to produce the GP's likelihood function input matrix
    # Signature: (<zfun_params>) -> Vector{Float64}
    zfun::Function
    # Parameters of the z function
    zfun_params::Array{Symbol}

    # The Gaussian Process object
    gp::AbstractGP
end


struct Parameter
    name::String
    longname::String
    default::Float64
    def_prior::UnivariateDistribution
    generate_link::Function
end

param_positive(p) = :(exp($p))
param_real(p) = :($p)

struct Likelihood
    z_inputs::Int
    params::Vector{Parameter}
    lik_expr::Expr
end

# Likelihood expressions have access to:
#  f, GP latent
#  z[k], likelihood input k
#  $param, likelihood parameter

data_likelihoods = Dict(
    "Gaussian" => Likelihood(
        0, # number of fields in z
        [Parameter("η", "StDev", 0.1, Exponential(1.0), param_positive))],
        η -> :(Normal(f, $η))
    ),
    "StudentT" => Likelihood(
        0,
        [Parameter("ν", "DOF", 4, Exponential(5.0), param_positive)),
         Parameter("σ", "Scale", 0.1, Exponential(1.0), param_positive))],
        ν, σ -> :(LocationScale(f, $σ, TDist($ν)))
    ),
    "Binomial" => Likelihood(
        1,
        [Parameter("u", "CompositionalMean", 0, Exponential(1.0), param_real)],
        u -> quote
            efu = exp(f + $u)
            Binomial(z[1], efu / (1.0 + efu))
        end
    ),
    "BetaBinomial" => Likelihood(
        1,
        [Parameter("u", "CompositionalMean", 0, Exponential(1.0), param_real),
         Parameter("r", "Overdispersion", 1, Exponential(1.0), param_positive))],
        u, r -> quote
            # Reparameterize as mean and variance (p, q)
            efu = exp(f + $u)
            p = efu / (1.0 + efu)
            # q = (z[i,1] - r) / (r - 1) # Target formula
            q = z[1] / $r # Approximation with no negatives/zero problems

            # Transform to (a, b) from (p, q)
            a = q * p
            b = q * (1 - p)

            BetaBinomial(z[1], a, b)
        end
    )
)

struct CovarianceFunction
    x_inputs::Int
    has_scalar::Bool
    params::Vector{Parameter}
    cf_expr::Function
end

# Covariance functions have access to:
#  i, j, sample numbers
#  $(t[k]), metadata index for metadata k
#  x[i,$(t[k])], parameter vector for sample i and metadata k
#  θ[$s+k], covariance function parameter k

covariance_functions = Dict(
    "Constant" => CovarianceFunction(
        0, true,
        [Parameter("σ2", "Variance", 1, Gamma(2, 0.5), i -> :(exp(θ[$i])))],
        (s, t) -> :(θ[$s+1])
    ),
    "Noise" => CovarianceFunction(
        0, false, [],
        (s, t) -> :(i == j ? 1.0 : 0.0)
    ),
    "Cat" => CovarianceFunction(
        1, false, [],
        (s, t) -> :(x[i,$(t[1])] == x[j,$(t[1])] ? 1.0 : 0.0)
    ),
    "OO" =>  CovarianceFunction(
        1, false,
        [Parameter("l", "LengthScale", 0.1, Exponential(1), i -> :(exp(θ[$i])))],
        (s, t) -> :(exp(-abs(x[i,$(t[1])] - x[j,$(t[1])]) / θ[$s+1]))
    ),
    "SExp" =>  CovarianceFunction(
        1, false, [],
        (s, t) -> quote
            dxs = (x[i,$(t[1])] - x[j,$(t[1])]) / θ[$s+1]
            exp(-0.5 * (dxs * dxs))
        end
    )
)
covariance_functions["1"] = covariance_functions["Constant"]


# Variable allocation for likelihood/covariance functions
struct VarAllocator
    var_idxs::Dict{String, Int}
    var_iscat::Dict{String, Bool}

    VarAllocator(varnames, variscat) = begin
        iscat = Dict{String, Bool}()
        for i in 1:length(varnames)
            iscat[varnames[i]] = variscat[i]
        end
        new(Dict{String, Int}(), iscat)
    end
end

function var_exists(alloc::VarAllocator, what::String)
    return what in alloc.var_iscat.keys()
end

function var_iscat(alloc::VarAllocator, what::String)
    return var_exists(alloc, what) && alloc.var_iscat[what]
end

function allocate(alloc::VarAllocator, what::String)
    # Allocate a variable and return its index
    if what in alloc.var_idxs.keys()
        return alloc.var_idxs[what]
    elseif var_exists(alloc, what)
        ix = length(alloc.var_idxs) + 1
        alloc.var_idxs[what] = ix
        return ix
    else
        error(@sprintf("Unknown variable: %s", what))
    end
end

function allocated(alloc::VarAllocator)
    # Get a list of allocated variables in the order they are allocated
    varlist = repeat("", inner=length(alloc.var_idxs))
    for (var, ix) in alloc.var_idxs
        varlist[ix] = var
    end
    return varlist
end

function generate_generators(alloc::VarAllocator)
    # Fun-named function which makes GP input generation functions from the
    # variables allocated in this allocator

    vars = allocated(alloc)

end


function sanitize_var_name(name)
    if !(isletter(name[1]) || name[1] == '_')
        name = string("X", name)
    end
    return replace(name, r"[^A-Za-z0-9_]" => "_")
end

function parse_gp_formula(formula::String, var_names::Array{String}, var_cat::Array{Bool})
    # Not as good as a real parser, but works for our purposes
    hasmatch, yformula, final, s = match_until_wnested(formula, Set([':', '~']))
    if !hasmatch
        error("Expected formula for y")
    end

    # Get potential variable names
    variablenames = [sanitize_var_name(string(name)) name in var_names]
    varname_set = Set{Symbol}(Symbol.(variablenames))

    # Evaluate y
    try
        yex = Meta.parse(yformula)
        yvars = Array{Symbol}(getvariables(yex))
        yfun = SyntaxTree.genfun(yex, yvars)
    catch err
        error(@sprintf("Error in Y formula: %s", string(err)))
    end

    # Parse the likelihood (optional)
    if final == ':'
        zalloc = VarAllocator(variablenames, var_cat)
        lik_expr, θ, θl_prior, θl_link_expr, θl_names, s = parse_lik(s, zalloc)
        s = chompw(s)
        if isempty(s) || s[1] != '~'
            error("Expected ~ after data likelihood.")
        end
        s = s[2:end]
    else
        lik = :(Normal(f, θ[1]))
        θl_prior = [data_likelihoods["Gaussian"].params[1].def_prior]
        θl_link = [data_likelihoods["Gaussian"].params[1].generate_link(1)]
    end
    if s[1] != '|'
        error("Expected ~|")
    end
    s = s[2:end]
    datalik = SyntaxTree.genfun(lik_expr, [:f, :z, :θ])
    θl_link = SyntaxTree.genfun(θl_link_expr, [:θ])

    # Parse the covariance function
    x = zeros(nrow(data), 0)
    θc = zeros(0)
    θc_prior = Array{UnivariateDistribution, 1}()
    θc_link = Array{Expr, 1}()

    cf_ex, s = parse_cf_expression(s, x, θc, θc_prior, θc_link, data, variablenames, variablesymbols, true)

    # Turn the covariance function into an actual Julia function
    cf = SyntaxTree.genfun(cf_ex, [:i, :j, :θ])

    return GP(dataloglik, (θl_prior...,), (θl_link...,),
              cf, (θc_prior...,), (θc_link...,))
end

function parse_lik(s::String, zalloc::VarAllocator)
    # lik_expr, θ, θl_prior, θl_link_expr, θl_names, s

    # Seek to first character
    os = s;
    s = chompw(s);
    if isempty(s) || !isletter(s[1])
        error("Expected likelihood")
    end

    # Get the likelihood name
    word, s = chomp_name(s)
    if !(name in keys(data_likelihoods))
        error(@sprintf("Unknown data likelihood: %s", name))
    end
    s = chompw(s)
    lik = data_likelihoods[name]

    # Set up parameters
    params = lik.params
    θ = [p.default for p in params]
    θ_prior = [p.def_prior for p in params]

    if isempty(s) || s[1] != '('
        if lik.z_input > 0
            error("Expected input metadata")
        elseif s[1] != '~'
            error("Expected ~")
        end
    else
        s = parse_params(s, zalloc, θ, θ_prior, lik.z_input)
    end

    # Remove fixed priors
    fixed = isa.(θ_prior, DeltaDistribution)
    θ_idx = cumsum(.!fixed)
    θ_param_exprs = [fixed[i] ? θ_prior[i].value : :(θ[$(θ_idx[i])]) for i in eachindex(fixed)]
    θ = θ[.!fixed]
    θ_prior = θ_prior[.!fixed]

    # Evaluate link and names
    θ_link_expr = Expr(:tuple, [p.link(i) for (i, p) in enumerate(params[.!fixed])]...)
    θ_names = [@sprintf("θl[%s]", p.name) for p in params[.!fixed]]

    # Generate the expression
    lik_expr = lik.lik_expr(θ_param_exprs...)

    return lik_expr, θ, θ_prior, θ_link_expr, θ_names, s
end

function parse_cf_expression(s, x, θ, θ_prior, θ_link, data, variablenames, variablesymbols, toplevel)

    cf_ex = :()
    overallneedsparam = true
    needsparam = true
    isprod = false

    while true
        s = chompw(s)
        if isempty(s)
            error("Expected: covariance function formula");
        elseif isletter(s[1]) || s[1] == '1'
            # Covariance function name or variable name
            if s[1] == '1'
                word = "1"
                s = s[2:end]
            else
                word, s = chomp_name(s)
            end

            if word in keys(covariance_functions)
                # This is a known covariance function
                cf = covariance_functions[word]

                # Set up default parameters and priors
                cf_θ = map(p -> p.default, cf.params)
                cf_θ_prior = map(p -> p.def_prior, cf.params)

                # Parse the parameters and metadata use
                if isempty(s) || s[1] != '('
                    if cf.x_inputs > 0
                        error("Expected input metadata")
                    end
                else
                    s, xx = parse_params(s, data, variablenames, variablesymbols, cf_θ, cf_θ_prior, cf.x_inputs)
                end
            elseif word in variablenames
                # Variable name
                idx = findfirst(word .== variablenames)
                xx = data[idx]

                if !isa(xx, Number)
                    xx = collapse_to_numeric(xx)
                    cf = covariance_functions["Cat"]
                else
                    cf = covariance_functions["Linear"]
                end

                # Use defaults for the cf
                cf_θ = map(p -> p.default, cf.params)
                cf_θ_prior = map(p -> p.def_prior, cf.params)
            else
                error(@sprintf("Unknown identifier in formula: %s", word))
            end

            # Merge new parameters used by the cf
            cf_s = length(θ)
            cf_θ_link = map(i -> cf.params[i].generate_link(length(θ_link) + i), eachindex(cf.params))
            θ = vcat(θ, cf_θ)
            θ_prior = vcat(θ_prior, cf_θ_prior)
            θ_link = vcat(θ_link, cf_θ_link)

            # Merge new metadata used by the cf into x
            cf_t = zeros(cf.x_inputs)
            for i = 1:cf.x_inputs
                hit = false
                for j = 1:ncol(x)
                    if all(x[:,j] .== xx[:,i])
                        cf_t[i] = j
                        hit = true
                        break
                    end
                end
                if !hit
                    x = hcat(x, xx[:,i])
                    cf_t[i] = ncol(x)
                end
            end

            # Create the covariance function expression
            new_ex = cf.cf_expr(cf_s, cf_t)
            needsparamsub = cf.has_scalar
        elseif s[1] == '('
            new_ex, s, needsparamsub = parse_cf_expression(
                s[2:end], x, θ, θ_prior, data, variablenames, variablesymbols)
            new_ex = :(($new_ex))
        else
            error(@sprintf("Unexpected symbol in covariance function formula: %s", s[1]))
        end

        needsparam = needsparam & needsparamsub
        overallneedsparam = overallneedsparam & !needsparamsub

        # Append the new expression
        if cf_ex == :()
            cf_ex = new_ex
        elseif isprod
            cf_ex = :($cf_ex * $new_ex)
        else
            cf_ex = :($cf_ex + $new_ex)
        end

        function add_magnitude()
            cf = covariance_functions["Constant"]
            cf_s = length(θ)
            θ = vcat(θ, map(p -> p.default, cf.params))
            θ_prior = vcat(θ_prior, map(p -> p.def_prior, cf.params))
            θ_link = vcat(θ_link, map(i -> cf.params[i].generate_link(length(θ_link) + i), eachindex(cf.params)))
            new_ex = cf.cf_expr(cf_s, [])
            cf_ex = :($cf_ex * $new_ex)
        end

        # Continue the expression
        s = chompw(s)
        if isempty(s) || s[1] == ')'
            if !overallneedsparam && needsparam || toplevel && needsparam
                # Add an explicit magnitude parameter to this term if this
                # subexpression contains sums
                add_magnitude()
                overallneedsparam = false
            end
            if !isempty(s)
                if toplevel
                    error("Unexpected )")
                end
                s = s[2:end]
            end
            break
        elseif s[1] == '+'
            if needsparam
                # This term needs an explicit magnitude parameter
                add_magnitude()
            end

            s = s[2:end]
            needsparam = true
            overallneedsparam = false
            isprod = false
        elseif s[1] == '*'
            s = s[2:end]
            isprod = true
        else
            error(@sprintf("Unexpected \'%s\' in covariance function formula", s[1]))
        end
    end

    return cf_ex, s, overallneedsparam
end

function parse_params(s, zalloc, θ, θ_prior, lik.z_input)
    xi = 1
    while xi <= n_inputs
        hasmatch, formula, final, s = match_until_wnested(s, Set([':', '~', ')']))

        if !hasmatch
            error("Expected )")
        end

        try
            xex = Meta.parse(formula)
            xfun = SyntaxTree.genfun(xex, variablesymbols)
            x[:,xi] = collapse_to_numeric(map(i -> xfun(map(j->data[xi,j], 1:ncol(data))...), 1:nrow(data)))
        catch err
            error(@sprintf("Error in metadata formula: %s", string(err)))
        end

        xi += 1
        if xi <= n_inputs
            if final != ','
                error("Expected ,")
            end
        elseif final == ')'
            s = final * s
        elseif final == ','
            error("Expected ; or )")
        end
    end

    if s[1] != ')'
        hasmatch, priors, final, s = match_until_wnested(s, Set([')']))

        if !hasmatch
            error("Expected )")
        end

        try
            θ_prior_over = eval(Meta.parse(priors))
        catch err
            error(@sprintf("Error in priors: %s", string(err)))
        end

        if length(θ_prior_over) > length(θ_prior)
            error("Too many parameters")
        end

        for i = eachindex(θ_prior_over)
            θ_prior[i] = θ_prior_over[i]
        end
    end

    return s, x
end

function collapse_to_numeric(x)
    if isa(x[1], Bool)
        return map(xi -> xi ? 1.0 : 0.0, x)
    elseif isa(x[1], Number)
        return x
    end

    # Assume x is a factor.. translate to levels
    # Could be more efficient...
    levels = unique(x)
    return map(xi -> findfirst(xi .== levels), x)
end

function chompw(s)
    i = 1
    while i <= length(s) && isspace(s[i])
        i += 1
    end
    return i > 1 ? s[i:end] : s
end

function chomp_name(s)
    i = 1
    while i <= length(s) && isletter(s[i])
        i += 1
    end
    return s[1:(i-1)], s[i:end]
end

function chomp_number(s)
    i = 1
    mt = match(r"^(\-?[0-9]*\.?[0-9]*[eE]?\-?[0-9]*)(.*)$")
    if isa(mt, Nothing)
        return parse(Float64, mt[1]), mt[2]
    end
    return nothing, s
end

function match_until_wnested(s, what::Set{Char})
    i = 1
    paren_count = 0
    while true
        if i > length(s)
            # Could not find what
            return false, "", "", s
        elseif s[i] == '('
            paren_count += 1
        elseif s[i] == ')' && paren_count > 0
            paren_count -= 1
        elseif paren_count == 0 && s[i] in what
            return true, s[1:(i-1)], s[i], s[(i+1):end]
        end
        i += 1
    end
end

# Get the set of variables used by an expression
getvariables(ex, vars::Set{Symbol}) = Set{Symbol}()
getvariables(ex::Symbol, vars::Set{Symbol}) = ex in vars ? Set{String}([ex]) : Set{String}()
getvariables(ex::Expr, vars::Set{Symbol}) = foldl(union, getvariables.(ex.args))

function gp_inputs(pf::ParsedGPFormula, data::DataFrame)
    # Utility function to get the GP input variables from a data frame and a
    # parsed GP formula

    # Evaluate x
    xdata = data[pf.xfun_params]
    x = vcat([pf.xfun(((xdata[i,j] for j in 1:size(xdata)[2])...,)) for i in 1:size(xdata)[1]]...)

    # Evaluate y
    ydata = data[pf.yfun_params]
    y = [pf.yfun(((ydata[i,j] for j in 1:size(ydata)[2])...,)) for i in 1:size(ydata)[1]]

    # Evaluate z
    zdata = data[pf.zfun_params]
    z = vcat([pf.zfun(((zdata[i,j] for j in 1:size(zdata)[2])...,)) for i in 1:size(zdata)[1]]...)

    return x, y, z
end
