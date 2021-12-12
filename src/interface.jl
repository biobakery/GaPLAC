struct Spec
    response::Symbol
    lik
    formula
end

likelihood(gps::Spec) = gps.lik
response(gps::Spec) = gps.response
formula(gps::Spec) = gps.formula


function gp_spec(formula::AbstractString)
    spl1 = findfirst(==(':'), formula)
    spl2 = last(findfirst(==('~'), formula))
    isnothing(spl2) && throw(ArgumentError("Invalid formula specification"))
    barind = nextind(formula, spl2)
    formula[barind] == '|' || throw(ArgumentError("Invalid formula specification"))

    if isnothing(spl1) || spl1 > spl2
        lik = GaPLAC.Gaussian()
        spl1 = prevind(formula, spl2)
    else
        spl1 < spl2 || throw(ArgumentError("Invalid formula specification"))
        lik = strip(formula[nextind(formula, spl1):prevind(formula, spl2)])
        lik = isempty(lik) ? GaPLAC.Gaussian() : GaPLAC.eval(Meta.parse(lik))
    end

    resp = strip(formula[1:prevind(formula, spl1)])
    
    gp = strip(formula[nextind(formula, barind):end])
    gp = GaPLAC.eval(Meta.parse(gp))

    return Spec(Symbol(resp), lik, gp)
end

function make_gp(spec::Spec)
    kern, vars = GaPLAC.kernel(GaPLAC.formula(spec))
    kernels = GaPLAC._walk_kernel(kern)
    length(vars) == length(kernels) || error("Something went wrong with equation parsing, number of variables should == number of kernels")
    return (GP(kern), vars)
end


function getatrange(at::AbstractString, vars=Symbol[])
    atdict = Dict{Symbol, Any}()
    ats = Meta.parse(at)
    @debug "Passed ranges" ats

    exprs = ats.head == :toplevel ? ats.args : [ats]
    
    for expr in exprs
        expr.head == Symbol("=") || error("Only assignments allowed in `--at` argument")
        var = expr.args[1]
        val = GaPLAC.eval(expr.args[2])
        atdict[var] = val
    end
    
    for (i, var) in enumerate(vars)
        if !in(var, keys(atdict))
            atdict[var] = GaPLAC._default_range(kernels[i])
        end
    end
    return atdict
end

@testset "Interface" begin
    @testset "Formula Parsing" begin
        spec1 = gp_spec("y ~| SqExp(:t)")
        @test likelihood(spec1) isa GaPLAC.Gaussian
        @test response(spec1) == :y
        @test formula(spec1) isa GaPLAC.GPCompnent
        @test formula(spec1) isa GaPLAC.SqExp
        
        spec2 = gp_spec("bug ~| SqExp(:t) + Linear(:x)")
        @test likelihood(spec2) isa GaPLAC.Gaussian
        @test response(spec2) == :bug
        @test formula(spec2) isa GaPLAC.GPCompnent
        @test formula(spec2) isa GaPLAC.GPOperation

        spec3 = gp_spec("bug ~| SqExp(:t) * Cat(:g) + Linear(:x)")
        @test likelihood(spec3) isa GaPLAC.Gaussian
        @test response(spec3) == :bug
        @test formula(spec3) isa GaPLAC.GPCompnent
        @test formula(spec3) isa GaPLAC.GPOperation
    end

    @testset "Ranges" begin
        at1 = getatrange("x = rand(Uniform(-5,5), 50)")
        @test haskey(at1, :x)
        @test length(at1[:x]) == 50
        low,high = extrema(at1[:x])
        @test -5 < low < high < 5
        
        at2 = getatrange("thing = rand(Normal(0,1), 100)")
        @test haskey(at2, :thing)
        @test length(at2[:thing]) == 100
        @test -0.5 < mean(at2[:thing]) < 0.5
    end
end