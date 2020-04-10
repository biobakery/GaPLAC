
abstract type AbstractChains
end

mutable struct Chains <: AbstractChains
	# Main dataframe
	df::DataFrame

	Chains() = new(DataFrame())
	Chains(d::DataFrame) = new(d)
end

# Maybe overload Base.getindex()?
function getrecords(c::Chains, what)
    # Get data for a sample
    c.df[!, what] # this returns a view, not a copy. Can use `c.df[:, what]` to get copy
end

function record!(c::Chains, what::Symbol, value::Float64)
	# Get indices of the new record
	ix = findfirst(isequal(what), names(c.df))

    # Add fields that we haven't seen yet
    if isa(ix, Nothing)
        # when the dataframe is totally empty, this adds a 1-row column w/NaN, is that intended?
        # otherwise the last line of this function throws an error when the DataFrame is empty
        nrow(c.df) == 0 ? c.df[!, what] = [NaN] : c.df[!,what] .= NaN
        ix = findfirst(isequal(what), names(c.df)) # since you added column, this could be `size(c.df, 2)` instead
    end

	# Fill in the values
	c.df[end,ix] = value
end

function thin(c::Chains, burnin::Integer, thinning::Integer)
	# Chain thinning
	# Does it matter that this thinnning isn't random?
	Chains(c.df[(burnin+1):thinning:end, :])
end

function newsample!(c::Chains)
	# Begin a new sample
	push!(c.df, repeat([NaN], size(c.df,2)))
end

# IO
function read_chains(filename; tsv::Bool = false)
	Chains(CSV.read(filename, delim=tsv ? '\t' : ','))
end
function write_chains(c::Chains, filename; tsv::Bool = false, append::Bool = false)
	CSV.write(filename, c.df, delim=tsv ? '\t' : ',', append=append)
end
