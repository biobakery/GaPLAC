struct DirectGP <: AbstractGP
	# The covariance function
	# Signature: (x_i::Vector, x_j::Vector, self::Boolean, θc::Vector) -> Real
	cf::Function # Covariance function used for prediction
	cftr::Function # Covariance function used for training
	# Covariance function hyperparameter transformation
	θc_link::Function
	# Prior on covariance function hyperparameters
	θc_prior::Vector{UnivariateDistribution}
	# Covariance function hyperparameter names
	θc_names::Vector{String}
	# Jitter
	jitter::Float64
end


function θ(gp::DirectGP, ϕ)
	return gp.θc_link(ϕ)
end

function ∇θ(gp::DirectGP, ϕ)
	θcd = θ(gp, [ForwardDiff.Dual(ϕi, 1.) for ϕi in ϕ])
	return ForwardDiff.value.(θcd), ForwardDiff.partials.(θcd,1)
end

function covariance_function!(C::Matrix, gp::DirectGP, θc, x)
	for i in 1:size(x,1)
		for j in i:size(x,1)
			if i==j
				C[i,i] = gp.cftr(x[i,:], x[j,:], true, θc) + gp.jitter
			else
				C[i,j] = C[j,i] = gp.cftr(x[i,:], x[j,:], false, θc)
			end
		end
	end
end

function covariance_function!(C::Matrix, gp::DirectGP, θc, x1, x2)
	for i in 1:size(x1,1)
		for j in 1:size(x2,1)
			C[i,j] = gp.cf(x1[i,:], x2[j,:], false, θc)
		end
	end
end

function record!(gp::DirectGP, chains::Chains, ϕ)
	# Transform parameters
	θc = θ(gp, ϕ)

	for i in eachindex(θc)
		record!(chains, Symbol(gp.θc_names[i]), θc[i])
	end
end

function invlink(gp::DirectGP, θ_target)
	# Bisecting search for ϕ which produces θ
	minv = repeat([-1000.0], inner=length(θ_target))
	maxv = repeat([1000.0], inner=length(θ_target))
	n = 0
	while any((maxv .- minv) .> 4.0 * max(eps.(minv), eps.(maxv))) && n < 64
		midv = (maxv .+ minv) ./ 2.0
		θ_mid = vcat(θ(gp, midv)...)
		gt = θ_mid .> θ_target
		maxv[gt] .= midv[gt]
		minv[.!gt] .= midv[.!gt]
		n += 1
	end

	return maxv
end

function unrecord(gp::DirectGP, chains::Chains, ix::Int)
	# Build the target θ vector
	θ_target = [getrecords(chains, Symbol(name))[ix] for name in gp.θc_names]

	return invlink(gp, θ_target)
end

function predict(gp::DirectGP, mcmc, x, y, z, x2, z2; quantiles)
	# Transform parameters
	θc = θ(gp, unrecord(gp, mcmc, size(mcmc.df,1)))

	# Evaluate covariance matrices
	xCx = zeros(size(x, 1), size(x, 1))
	covariance_function!(xCx, gp, θc, x)
	xCy = zeros(size(x2, 1), size(x, 1))
	covariance_function!(xCy, gp, θc, x2, x)
	yCy = zeros(size(x2, 1), size(x2, 1))
	covariance_function!(yCy, gp, θc, x2)

	# Conditional distribution
	μf = (xCy * (xCx \ y))[:,1]
	σ2f = diag(yCy) - diag(xCy * (xCx \ transpose(xCy)))
	
	Q = zeros(size(x2, 1), length(quantiles))
	for i in 1:length(μf)
		Q[i,:] = quantile.(Normal(μf[i], sqrt(σ2f[i])), quantiles)
	end

	return μf, σ2f, Q, μf, σ2f
end

function samplegp(gp::DirectGP, ϕ, x, y, z, x2, z2)
	# Transform parameters
	θc = θ(gp, ϕ)

	# Evaluate covariance matrices
	xCx = zeros(size(x, 1), size(x, 1))
	covariance_function!(xCx, gp, θc, x)
	xCy = zeros(size(x2, 1), size(x, 1))
	covariance_function!(xCy, gp, θc, x2, x)
	yCy = zeros(size(x2, 1), size(x2, 1))
	covariance_function!(yCy, gp, θc, x2)

	# Conditional distribution
	μf2 = xCy * (xCx \ y)
	Σf2 = yCy - xCy * (xCx \ transpose(xCy))

	# Sample the gp at the target points
	# Uses eigendecomposition to allow positive semi-definite covariance
	# matrices (could be optimized)
	E = eigen(Σf2)
	L = real.(E.vectors) * diagm(0 => sqrt.(max.(real.(E.values), 0.0)))
	y_samp = L * randn(length(μf2)) .+ μf2

	return y_samp, y_samp
end

function logposterior(gp::DirectGP, ϕ, x, y, z)
	# Separate and transform parameters
	θc, dθc_dϕc = ∇θ(gp, ϕ)

	# Evaluate the covariance function at θ
	C = zeros(size(x, 1), size(x, 1))
	covariance_function!(C, gp, θc, x)

	# Evaluate loglik
	ll = logpdf(MvNormal(zeros(size(x,1)), C), y)

	# Log prior
	lp = sum(logpdf.(gp.θc_prior, θc) .+ log.(dθc_dϕc))

	# Return log prior and log lik separately
	return lp, ll
end
