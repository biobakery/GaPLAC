
struct AbstractMCMC
end

mutable struct MCMCStats
	# Total time taken by the sampler
    total_time_s::Float64
	# Accepted samples
	accepts::Int
	# Rejected samples
    rejects::Int
	# Total function evaluations
    evaluations::Int
	# Total gradient function evaluations
    grad_evaluations::Int

    new = MCMCStats(0.0, 0, 0, 0, 0)
end

struct UnivariateSliceMCMC :< AbstractMCMC
	# Target log density
	# Signature: (x::Vector{Float64}) -> Float64
	f::Function
	# Parameter widths
	w::Vector{Float64}
end


function sample!(x::AbstractVector, mcmc::UnivariateSliceMCMC, stats::MCMCStats)
    # TODO: Adapt code
end

# OLD CODE:
function uvslicesamp!(x, w, llfun::Function; N=1)
	ll_x0 = llfun(x)
	nreject = 0

	for i = 1:N
		for k = eachindex(x)
			# Set up the window
			ub = w
			lb = -w
			sll = ll_x0 + log(rand())

			# Shrink until we're inside the slice
			dx = 0.0 .* x
			lln = 0.0
			for j = 1:20
				# Pick a point inside the window
				dx[k] = lb + rand() * (ub - lb)

				# Are we inside?
				lln = llfun(x .+ dx)
				lln >= sll && break # Yup

				# Nope - shrink
				nreject += 1
				if dx[k] < 0
					lb = dx[k]
				else
					ub = dx[k]
				end
			end
			if lln < sll
				# Failed
				continue
			end

			# Move
			x .= x .+ dx
			ll_x0 = lln
		end
	end

	return ll_x0, nreject
end
