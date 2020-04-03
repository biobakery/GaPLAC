
abstract type AbstractMCMC end

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

	MCMCStats() = new(0.0, 0, 0, 0, 0)
end

struct UnivariateSliceMCMC <: AbstractMCMC
	# Target log density
	# Signature: (x::Vector{Float64}) -> Float64
	f::Function
	# Parameter widths
	w::Vector{Float64}
	# Number of samples to generate
	N::Int
end


function mcmcsample!(x::AbstractVector, mcmc::UnivariateSliceMCMC, stats::MCMCStats = MCMCStats())
	start_t = time_ns()
	y_x0 = mcmc.f(x)
	stats.evaluations += 1

	for i = 1:mcmc.N
		for k = eachindex(x)
			# Set up the window
			lb, ub = -mcmc.w[k], mcmc.w[k]
			sy = y_x0 + log(rand())

			# Shrink until we're inside the slice
			dx = zeros(length(x))
			yn = 0.0
			for j = 1:20
				# Pick a point inside the window
				dx[k] = lb + rand() * (ub - lb)

				# Are we inside?
				yn = mcmc.f(x .+ dx)
				stats.evaluations += 1
				yn >= sy && break # Yup

				# Nope - shrink
				stats.rejects += 1
				if dx[k] < 0
					lb = dx[k]
				else
					ub = dx[k]
				end
			end
			if yn < sy
				# Failed
				continue
			end

			# Move
			stats.accepts += 1
			x .= x .+ dx
			y_x0 = yn
		end
	end

	stats.total_time_s += (time_ns() - start_t) / 1.e9
	return y_x0
end
