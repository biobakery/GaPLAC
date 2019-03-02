
# round(y * filtered_reads) : BetaBinomialR(filtered_reads) ~|
#   Cat(subject) * OU(time; l=Gamma(1)) +
#   Linear(diet)

# TODO: Errors from these functions could provide better context so that
# finding the actual problem in the formula is easier

using DataFrames
using Distributions
using Printf

struct Delta{T<:Real} <: ContinuousUnivariateDistribution
    # The single value of the distribution
    value::T
end

struct GPComponent
    # Describes a part of a GP covariance function

    # The indices of the x variables used by the component
    x_varindices::Array{Int,1}
    # The covariance function for the component
    cf::Function
    cf_expr::Any
end

struct ParsedGPFormula
    # Function to produce the GP's covariance function input matrix
    # Signature: (<xfun_params>) -> Vector{Float64}
    xfun::Function
    # Parameters of the x function
    xfun_params::Array{Symbol,1}
    # Variables used in each x variable
    xvars::Array{Set{Symbol},1}
    # Names of the x matrix columns
    xnames::Array{String,1}

    # Function to produce the GP's output
    # Signature: (<yfun_params>) -> Number
    yfun::Function
    # Parameters of the y function
    yfun_params::Array{Symbol}
    # Name of the y values
    yname::String

    # Function to produce the GP's likelihood function input matrix
    # Signature: (<zfun_params>) -> Vector{Float64}
    zfun::Function
    # Parameters of the z function
    zfun_params::Array{Symbol,1}
    # Variables used in each x variable
    zvars::Array{Set{Symbol},1}
    # Names of the z matrix columns
    znames::Array{String,1}

    # The Gaussian Process object
    gp::AbstractGP
    # Components of the GP
    components::Array{GPComponent,1}

    # Starting parameters
    θc::Array{Float64,1}
    θl::Array{Float64,1}
end

function gp_predcomponent(formula::ParsedGPFormula, comp::Int)
    # Returns a GP which will produce predictions for a given component
    gp = copy(formula.gp)
    gp.cf = formula.components[comp].cf

    return gp
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
    lik_expr::Function
end

# Likelihood expressions have access to:
#  f, GP latent
#  z[k], likelihood input k
#  $param, likelihood parameter

data_likelihoods = Dict(
    "None" => Likelihood(0, [], () -> :(Normal(f, 1e-6))),
    "Gaussian" => Likelihood(
        0, # number of fields in z
        [Parameter("η", "StDev", 0.1, Exponential(1.0), param_positive)],
        η -> :(Normal(f, $η))
    ),
    "StudentT" => Likelihood(
        0,
        [Parameter("ν", "DOF", 4, Exponential(5.0), param_positive),
         Parameter("σ", "Scale", 0.1, Exponential(1.0), param_positive)],
        (ν, σ) -> :(LocationScale(f, $σ, TDist($ν)))
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
         Parameter("r", "Overdispersion", 1, Exponential(1.0), param_positive)],
        (u, r) -> quote
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
#  $xi, $xj, covariance input
#  same, whether the two samples are the same sample
#  $param, covariance function parameter
# Signature: (xi..., xj..., params...) -> Expr

covariance_functions = Dict(
    "Constant" => CovarianceFunction(
        0, true,
        [Parameter("σ2", "Variance", 1, Gamma(2, 0.5), param_positive)],
        (σ2) -> :($σ2)
    ),
    "Noise" => CovarianceFunction(
        0, false, [],
        () -> :(same ? 1.0 : 0.0)
    ),
    "Cat" => CovarianceFunction(
        1, false, [],
        (cati, catj) -> :($cati == $catj ? 1.0 : 0.0)
    ),
    "OU" =>  CovarianceFunction(
        1, false,
        [Parameter("l", "LengthScale", 0.1, Exponential(1), param_positive)],
        (ti, tj, l) -> :(exp(-abs($ti - $tj) / $l))
    ),
    "SExp" =>  CovarianceFunction(
        1, false,
        [Parameter("l", "LengthScale", 0.1, Exponential(1), param_positive)],
        (ti, tj, l) -> quote
            dt = ($ti - $tj) / $l
            exp(-(dt * dt))
        end
    ),
    "Periodic" =>  CovarianceFunction(
        1, false,
        [Parameter("l", "LengthScale", 0.1, Exponential(1), param_positive),
         Parameter("ϕ", "Period", 1, Exponential(1), param_positive)],
        (ti, tj, l, ϕ) -> quote
            cdt = cos(($ti - $tj) * (π / $ϕ))
            exp(-(cdt*cdt)/$l)
        end
    )
)


# Variable allocation for likelihood/covariance functions
struct VarAllocator
    var_idxs::Dict{Symbol, Int}
    var_iscat::Dict{Symbol, Bool}
    varset::Set{Symbol}

    VarAllocator(varnames, variscat) = begin
        iscat = Dict{Symbol, Bool}()
        for i in 1:length(varnames)
            iscat[Symbol(varnames[i])] = variscat[i]
        end
        new(Dict{Symbol, Int}(), iscat, Set(Symbol.(varnames)))
    end
end

var_exists(alloc::VarAllocator, what) = Symbol(what) in alloc.varset
var_iscat(alloc::VarAllocator, what) =
    var_exists(alloc, what) && alloc.var_iscat[Symbol(what)]

function allocate!(alloc::VarAllocator, what)
    # Allocate a variable and return its index
    what = Symbol(what)
    if what in keys(alloc.var_idxs)
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
    varlist = repeat([:+], inner=length(alloc.var_idxs))
    for (var, ix) in alloc.var_idxs
        varlist[ix] = var
    end
    return varlist
end

genfun(ex, args) = eval(Expr(:->, Expr(:tuple, args...), ex))
generate_generator(alloc::VarAllocator, ex::Expr) =
    genfun(ex, allocated(alloc))
generate_generator(alloc::VarAllocator, ex::Array{Any}) =
    generate_generator(alloc, Expr(:vect, ex...))


function sanitize_var_name(name)
    if !(isletter(name[1]) || name[1] == '_')
        name = string("X", name)
    end
    return replace(name, r"[^A-Za-z0-9_]" => "_")
end

function parse_gp_formula(formula::String, var_names::Array{String}, var_cat::Array{Bool})
    # Main entrypoint for parsing a GP formula
    # Uses a simple hand-rolled recursive descent parser

    # Not as good as a real parser, but works for our purposes
    hasmatch, yformula, final, s = match_until_wnested(formula, Set([':', '~']))
    if !hasmatch
        error("Expected formula for y")
    end

    # Get potential variable names
    # TODO: Move this outside of this function so that var names are universal
    variablenames = [sanitize_var_name(string(name)) for name in var_names]
    varname_set = Set(Symbol.(variablenames))

    # Evaluate y
    yfun = identity
    yvars = []
    yex = :()
    try
        yex = Meta.parse(yformula)
        yvars = [getvariables(yex, varname_set)...]
        yfun = genfun(yex, yvars)
    catch err
        error(@sprintf("Error in Y formula: %s", string(err)))
    end

    # Parse the likelihood (optional)
    zalloc = VarAllocator(variablenames, var_cat)
    if final == ':'
        lik_ex, zex, znames, θl, θl_prior, θl_link_ex, θl_names, s = parse_lik(s, zalloc)
        s = chompw(s)
        if isempty(s) || s[1] != '~'
            error("Expected ~ after data likelihood.")
        end
        s = s[2:end]
    else
        lik_ex, zex, znames, θl, θl_prior, θl_link_ex, θl_names = parse_lik("Gaussian", zalloc)
    end
    if s[1] != '|'
        error("Expected ~|")
    end
    s = s[2:end]

    # Parse the covariance function
    xalloc = VarAllocator(variablenames, var_cat)
    θc, θc_prior, θc_names, θc_link_ex =
        Array{Float64,1}(), Array{ContinuousUnivariateDistribution,1}(),
        Array{String,1}(), Array{Any,1}()
    xex, xnames = Array{Any,1}(), Array{String,1}()
    if occursin(r"^\s*0\s*$", s)
        cf_ex = :(0.)
        s = ""
        comps = Array{Component,1}()
    else
        cf_ex, s, needsparam, comps = parse_cf_expression(s, θc, θc_prior, θc_names, θc_link_ex, xex, xnames, xalloc, true)
    end

    # Turn the data likelihood into an actual Julia function
    datalik = genfun(lik_ex, [:f, :z, :θ])
    θl_link = genfun(Expr(:vect, θl_link_ex...), [:θ])
    zfun = generate_generator(zalloc, zex)

    # Turn the covariance function into an actual Julia function
    cf = genfun(cf_ex, [:x1, :x2, :same, :θ])
    θc_link = genfun(Expr(:vect, θc_link_ex...), [:θ])
    xfun = generate_generator(xalloc, xex)

    # Record parse results
    @info "GP formula interpretation" Observation=yex Lik_Inputs=znames Lik_Parameters=θl_names Lik_Start=θl Likelihood=lik_ex CF_Inputs=xnames CF_Parameters=θc_names CF_Start=θc Covariance=cf_ex

    # Generate the GP object
    gp = LaplaceGP(
        cf, cf, # Training and prediction cfs are the same
        θc_link, θc_prior, θc_names, 1e-9,
        datalik, θl_link, θl_prior, θl_names)

    return ParsedGPFormula(
        xfun, allocated(xalloc), [getvariables(ex, xalloc) for ex in xex], xnames,
        yfun, yvars, yformula,
        zfun, allocated(zalloc), [getvariables(ex, zalloc) for ex in zex], znames,
        gp, comps,
        θc, θl)
end

function gp_inputs(pf::ParsedGPFormula, data::DataFrame)
    # Utility function to get the GP input variables from a data frame and a
    # parsed GP formula

    if !isempty(data)
        n = size(data)[1]

        # Evaluate x
        x = x = zeros(n, 0)
        if !isempty(pf.xnames)
            xdata = data[pf.xfun_params]
            x = vcat([transpose(pf.xfun((xdata[i,j] for j in 1:size(xdata)[2])...)) for i in 1:size(xdata)[1]]...)
        end

        # Evaluate y if all fields are available
        y = []
        if all([(p in names(data)) for p in pf.yfun_params])
            ydata = data[pf.yfun_params]
            y = [pf.yfun(((ydata[i,j] for j in 1:size(ydata)[2])...,)) for i in 1:size(ydata)[1]]
        end

        # Evaluate z
        z = zeros(n, 0)
        if !isempty(pf.znames)
            zdata = data[pf.zfun_params]
            z = vcat([transpose(pf.zfun((zdata[i,j] for j in 1:size(zdata)[2])...)) for i in 1:size(zdata)[1]]...)
        end

        return x, z, y
    else
        # Separate codepath for empty data so no problems occur if data has no columns
        x = zeros(0, length(pf.xnames))
        y = zeros(0)
        z = zeros(0, length(pf.znames))
        return x, z, y
    end
end

function parse_lik(s::String, zalloc::VarAllocator)
    # lik_expr, zex, znames, θ, θl_prior, θl_link_expr, θl_names, s

    # Seek to first character
    os = s;
    s = chompw(s);
    if isempty(s) || !isletter(s[1])
        error("Expected likelihood")
    end

    # Get the likelihood name
    name, s = chomp_name(s)
    if !(name in keys(data_likelihoods))
        error(@sprintf("Unknown data likelihood: %s", name))
    end
    s = chompw(s)
    lik = data_likelihoods[name]

    # Set up parameters
    params = lik.params
    θ = [p.default for p in params]
    θ_prior = convert(Array{ContinuousUnivariateDistribution}, [p.def_prior for p in params])

    if isempty(s) || s[1] == '~'
        if lik.z_inputs > 0
            error("Expected input metadata")
        end
        zex = []
        znames = []
    elseif s[1] == '('
        zex, znames, s = parse_params(s[2:end], zalloc, params, θ, θ_prior, lik.z_inputs)
    else
        error("Expected ( or ~")
    end

    # Remove fixed priors
    fixed = isa.(θ_prior, Delta)
    θ_idx = cumsum(.!fixed)
    θ_param_exprs = [fixed[i] ? θ_prior[i].value : :(θ[$(θ_idx[i])]) for i in eachindex(fixed)]
    θ = θ[.!fixed]
    θ_prior = θ_prior[.!fixed]

    # Evaluate link and names
    θ_link_expr = [p.generate_link(:(θ[$i])) for (i, p) in enumerate(params[.!fixed])]
    θ_names = [@sprintf("θl[%s]", p.name) for p in params[.!fixed]]

    # Generate the expression
    lik_expr = lik.lik_expr(θ_param_exprs...)

    return lik_expr, zex, znames, θ, θ_prior, θ_link_expr, θ_names, s
end

function parse_cf_expression(s, θ, θ_prior, θ_names, θ_link_ex, xex, xnames, xalloc, toplevel)
    # Returns cf_ex, s, overallneedsparam
    # Modifies θ, θ_prior, θ_names, θ_link_ex, xex, xnames, xalloc

    cf_ex = :()
    cfprod_ex = :()
    overallneedsparam = true
    needsparam = true
    isprod = false
    cf_θ = []
    cf_θ_prior = []

    cf_xnames = Set{Symbol}()
    comp_xidx = Array{Int,1}()
    comps = Array{GPComponent,1}()

    while true
        s = chompw(s)
        if isempty(s)
            error("Expected: covariance function formula")
        elseif isletter(s[1]) || s[1] == '1'
            # Covariance function name or variable name
            if s[1] == '1'
                word = "Constant"
                s = s[2:end]
            else
                word, s = chomp_name(s)
            end

            if word in keys(covariance_functions)
                # This is a known covariance function
                cf = covariance_functions[word]

                # Set up default parameters and priors
                cf_θ = [p.default for p in cf.params]
                cf_θ_prior = convert(Array{ContinuousUnivariateDistribution,1}, [p.def_prior for p in cf.params])

                # Parse the parameters and metadata use
                if isempty(s) || s[1] != '('
                    if cf.x_inputs > 0
                        error("Expected input metadata")
                    end
                    cf_xex = []
                    cf_xnames = []
                else
                    cf_xex, cf_xnames, s = parse_params(s[2:end], xalloc, cf.params, cf_θ, cf_θ_prior, cf.x_inputs)
                end
            elseif var_exists(xalloc, word)
                # Variable name
                s = @sprintf("%s(%s)%s", var_iscat(xalloc, word) ? "Cat" : "Linear", word, s)
                continue
            else
                error(@sprintf("Unknown identifier in formula: %s", word))
            end

            # Remove fixed priors and merge parameters
            fixed = isa.(cf_θ_prior, Delta)
            θ_idx = length(θ_prior) .+ cumsum(.!fixed)
            push!(θ, cf_θ[.!fixed]...)
            push!(θ_prior, cf_θ_prior[.!fixed]...)
            push!(θ_link_ex, [p.generate_link(:(θ[$(θ_idx[i])])) for (i, p) in enumerate(cf.params[.!fixed])]...)
            push!(θ_names, [begin
                name = @sprintf("θc[%s]", p.name)
                i = 1
                while name in θ_names
                    i += 1
                    name = @sprintf("θc[%s_%d]", p.name, i)
                end
                name
            end for p in cf.params[.!fixed]]...)

            # Merge x variables
            push!(comp_xidx, (length(xex) .+ collect(1:length(cf_xex)))...)
            xbase = length(xex)
            push!(xex, cf_xex...)
            push!(xnames, cf_xnames...)

            # Build parameters for the cf expression
            cf_θ_param_exprs = [fixed[i] ? cf_θ_prior[i].value : :(θ[$(θ_idx[i])]) for i in eachindex(fixed)]
            x1_params = [:(x1[$(xbase + i)]) for i in 1:cf.x_inputs]
            x2_params = [:(x2[$(xbase + i)]) for i in 1:cf.x_inputs]

            # Build the cf expression
            new_ex = cf.cf_expr(x1_params..., x2_params..., cf_θ_param_exprs...)
            needsparamsub = !cf.has_scalar
        elseif s[1] == '('
            # Read the subexpression
            new_ex, s, needsparamsub = parse_cf_expression(
                s[2:end], θ, θ_prior, θ_names, θ_link_ex, xex, xnames, xvars, xalloc, false)
        else
            error(@sprintf("Unexpected symbol in covariance function formula: %s", s[1]))
        end

        # Keep track of whether we need to add a scalar parameter to this component
        needsparam = needsparam & needsparamsub
        overallneedsparam = overallneedsparam & !needsparamsub

        # Aggregate products
        cfprod_ex = isprod ? :($cfprod_ex * $new_ex) : new_ex
        s = chompw(s)
        if !isempty(s) && s[1] == '*'
            s = s[2:end]
            isprod = true
            continue
        end
        isprod = false

        # Next term is not a continuation of the product, so make sure we add
        # a magnitude term if we need to
        if (isempty(s) || s[1] == ')') && (!overallneedsparam && needsparam || toplevel && needsparam) || needsparam
            s = "1 " * s
            isprod = true
            continue
        end

        # Track components
        toplevel && push!(comps, GPComponent(comp_xidx,
            genfun(cfprod_ex, [:x1, :x2, :same, :θ]), cfprod_ex))

        # Aggregate sums
        cf_ex = cf_ex == :() ? cfprod_ex : :($cf_ex + $cfprod_ex)

        # Continue the expression
        if isempty(s) || s[1] == ')'
            if !isempty(s)
                if toplevel
                    error("Unexpected )")
                end
                s = s[2:end]
            end
            break
        elseif s[1] == '+'
            s = s[2:end]
            needsparam = true
        else
            error(@sprintf("Unexpected \'%s\' in covariance function formula", s[1]))
        end
    end

    return cf_ex, s, overallneedsparam, comps
end

function parse_params(s, varalloc, θ_params, θ, θ_prior, nvars)
    # Parses a parameter set: (variables; parameter values/priors)
    # Returns var_ex, var_names, s
    # Modifies θ, θ_prior

    xi = 1
    var_ex = convert(Array{Any,1}, repeat([:()], inner=nvars))
    var_names = repeat([""], inner=nvars)
    while xi <= nvars
        # Read variables
        hasmatch, formula, final, s = match_until_wnested(s, Set([',', ';', ')']))

        if !hasmatch
            error("Expected )")
        end

        try
            var_names[xi] = formula
            var_ex[xi] = Meta.parse(formula)
            vars = getvariables(var_ex[xi], varalloc)
            varix = [allocate!(varalloc, var) for var in vars]
        catch err
            error(@sprintf("Error in metadata formula: %s", string(err)))
        end

        xi += 1
        if xi <= nvars
            if final != ','
                error("Expected ,")
            end
        elseif final == ')'
            return var_ex, var_names, s
        elseif final == ','
            error("Expected ; or )")
        end
    end

    # Hyperparameters
    hasmatch, priors, final, s = match_until_wnested(s, Set([')']))

    if !hasmatch
        error("Expected )")
    end

    unnamed_i = 1
    ss = chompw(priors * ")")
    while ss != ""
        hasmatch, prior, final, ss = match_until_wnested(ss, Set([',', ')']))

        if !hasmatch
            # This should be impossible
            error("Expected , or )")
        end

        m = match(r"^(?<name>[^=]+?) *= *(?<expr>.*)$", prior)
        if isa(m, Nothing)
            # Whole thing is an expression
            if unnamed_i > length(θ_params)
                error("Too many parameters")
            end
            prior_idx = unnamed_i
            expr = Meta.parse(prior) # TODO: Errors
            unnamed_i += 1
        else
            # name = value
            prior_idx = findfirst([p.name == m[:name] for p in θ_params])
            if isa(prior_idx, Nothing)
                prior_idx = findfirst([p.longname == m[:name] for p in θ_params])
            end
            if isa(prior_idx, Nothing)
                error(@sprintf("Unknown parameter: %s", m[:name]))
            end
            expr = Meta.parse(m[:expr]) # TODO: Errors
        end

        value = eval(expr)
        if isa(value, Number)
            # This is a single number - set the value and set it as a Delta
            θ[prior_idx] = value
            θ_prior[prior_idx] = Delta(value)
        elseif isa(value, UnivariateDistribution)
            # Set the prior to this distribution
            θ_prior[prior_idx] = value
            if isa(value, Delta)
                θ[prior_idx] = value
            end
        else
            error(@sprintf("Invalid value for %s: %s", θ_params[prior_idx].longname, string(typeof(value))))
        end
    end

    return var_ex, var_names, s
end

# Some utility functions

function chompw(s)
    # Chomp whitespace
    i = 1
    while i <= length(s) && isspace(s[i])
        i += 1
    end
    return i > 1 ? s[i:end] : s
end

function chomp_name(s)
    # Chomp a name
    i = 1
    while i <= length(s) && isletter(s[i])
        i += 1
    end
    return s[1:(i-1)], s[i:end]
end

function chomp_number(s)
    # Chome a number
    i = 1
    mt = match(r"^(\-?([0-9]+\.?[0-9]*|\.[0-9]+)([eE]\-?[0-9]+)?)(.*)$", s)
    if !isa(mt, Nothing)
        return parse(Float64, mt[1]), mt[4]
    end
    return nothing, s
end

function match_until_wnested(s, what::Set{Char})
    # Match a string until a terminal character is reached, balancing paired parens
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
getvariables(ex::Symbol, vars::Set{Symbol}) = ex in vars ? Set{Symbol}([ex]) : Set{Symbol}()
getvariables(ex::Expr, vars::Set{Symbol}) = foldl(union, [getvariables(exa, vars) for exa in ex.args])
getvariables(ex, vars::VarAllocator) = getvariables(ex, vars.varset)
