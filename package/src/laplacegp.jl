
using LinearAlgebra
using ForwardDiff
using Distributions
using FastGaussQuadrature

struct LaplaceGP <: AbstractGP
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

	# Data Likelihood
	# Signature: (f_i::Float64, z_i::Vector, θl::Vector) -> UnivariateDistribution
	datalik::Function
	# Data likelihood hyperparameter transformation
	θl_link::Function
	# Prior on data likelihood hyperparameters
	θl_prior::Vector{UnivariateDistribution}
	# Data likelihood hyperparameter names
	θl_names::Vector{String}
end


function θ(gp::LaplaceGP, ϕ)
	return gp.θl_link(ϕ[1:length(gp.θl_prior)]),
		gp.θc_link(ϕ[length(gp.θl_prior) .+ (1:length(gp.θc_prior))])
end

function ∇θ(gp::LaplaceGP, ϕ)
	θld, θcd = θ(gp, [ForwardDiff.Dual(ϕi, 1.) for ϕi in ϕ])
	return ForwardDiff.value.(θld), ForwardDiff.value.(θcd),
		ForwardDiff.partials.(θld,1), ForwardDiff.partials.(θcd,1)
end

function covariance_function!(C::Matrix, gp::LaplaceGP, θc, x)
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

function covariance_function!(C::Matrix, gp::LaplaceGP, θc, x1, x2)
	for i in 1:size(x1,1)
		for j in 1:size(x2,1)
			C[i,j] = gp.cf(x1[i,:], x2[j,:], false, θc)
		end
	end
end

function ∇2ll_f!(gp::LaplaceGP, ∇ll, ∇2ll, f, y, z, θl)
	# Log data likelihood and its first and second derivatives wrt f

	ll = 0.0
	for i = 1:length(f)
		@inbounds ydist = gp.datalik(ForwardDiff.Dual(ForwardDiff.Dual(f[i], 1.0), ForwardDiff.Dual(1.0, 0.0)), z[i,:], θl)
		v = logpdf(ydist, y[i])
		ll += ForwardDiff.value(ForwardDiff.value(v))
		@inbounds ∇ll[i] = ForwardDiff.partials(ForwardDiff.value(v), 1)
		@inbounds ∇2ll[i] = ForwardDiff.partials(ForwardDiff.partials(v, 1), 1)
	end

	return ll
end

function laplace_approx(gp::LaplaceGP, L::LowerTriangular, y, z, θl;
	lπtol = 0.01, maxN = 10)

	# Saddle-free Newton's method to find posterior mode for whitened latents fw
	f, fw = zeros(length(y)), zeros(length(y))
	∇ll_f, ∇lπ_fw = zeros(length(y)), zeros(length(y))
	∇2ll_f, ∇2ll_fw = zeros(length(y)), zeros(length(y), length(y))
	∇2ll_fw_temp = zeros(length(y), length(y))
	∇2lπ_fw, ∇2lπ_fw_reg = zeros(length(y), length(y)), zeros(length(y), length(y))
	dfw = zeros(length(y))
	α = 1.0
	lπ = -Inf

	for i in 1:maxN
		# Unwhiten
		mul!(f, L, fw)

		# Evaluate first and second derivatives of the likelihood at f
		ll = ∇2ll_f!(gp, ∇ll_f, ∇2ll_f, f, y, z, θl)
		@debug "W" ∇ll_f ∇2ll_f

		# Re-whiten Jacobian and Hessian
		mul!(∇lπ_fw, transpose(L), ∇ll_f)
		# Slow formula: ∇2ll_fw = transpose(L) * Diagonal(∇2ll_f) * L
		mul!(∇2ll_fw_temp, diagm(0 => ∇2ll_f), L)
		mul!(∇2lπ_fw, transpose(L), ∇2ll_fw_temp)

		# White noise prior on whitened latents
		lπ2 = ll - sum(fw.^2) / 2.0
		∇lπ_fw .-= fw
		for i in eachindex(fw)
			∇2lπ_fw[i,i] -= 1.0
		end

		# Slow if we're actually getting worse
		α = lπ2 < lπ ? 0.5 : 1.0
		lπ = lπ2
		@debug (i, lπ)

		E = eigen(∇2lπ_fw)
		Evals = E.values
		Evecs = E.vectors
		@debug "Eigs", Evals
		if !all(isreal.(Evals))
			# Sometimes, floating point accuracy makes the matrix slightly
			# assymetric, leading to imaginary eigenvalues. Discard them here.
			Evals = real.(Evals)
			Evecs = real.(Evecs)
		end
		if any(Evals .> -1.0)
			# We're in a part of the parameter space that has negative,
			# or where there is dangerously low curvature.
			# Apply a softabs transform here on eigenvalues to avoid saddle
			# points/minima and to regularize small values in the Hessian.
			λ = Evals .* coth.(Evals) # softabs with minimum value 1
			∇2lπ_fw_reg .= Evecs * diagm(0 => -λ) * transpose(Evecs)
			@warn "Applying softabs" λ
		else
			# We're in a safely-positive definite part of the parameter space
			∇2lπ_fw_reg .= ∇2lπ_fw

			# How much of a gain do we expect?
			#ldiv!(dfw, ∇2ll_fw, ∇ll_fw)
			dfw .= ∇2lπ_fw \ ∇lπ_fw
			@debug "dfw" ∇2lπ_fw ∇lπ_fw dfw
			lπgain = -transpose(dfw) * ∇2lπ_fw * dfw
			@debug "Gain" lπgain
			lπgain < lπtol && break # stop when we gain too little from another step
		end

		# Newton step
		if i < maxN
			fw .-= α * (∇2lπ_fw_reg \ ∇lπ_fw)
		else
			@warn "maxN reached in Newton steps"
		end
	end

	return fw, ∇2lπ_fw, lπ
end

function record!(gp::LaplaceGP, chains::Chains, ϕ)
	# Separate and transform parameters
	θl, θc = θ(gp, ϕ)

	for i in eachindex(θl)
		record!(chains, Symbol(gp.θl_names[i]), θl[i])
	end
	for i in eachindex(θc)
		record!(chains, Symbol(gp.θc_names[i]), θc[i])
	end
end

function invlink(gp::LaplaceGP, θ_target)
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

function extractθ(gp::LaplaceGP, chains::Chains, ix::Int)
	return [getrecords(chains, Symbol(name))[ix] for name in gp.θl_names],
		[getrecords(chains, Symbol(name))[ix] for name in gp.θc_names]
end

function unrecord(gp::LaplaceGP, chains::Chains, ix::Int)
	return invlink(gp, vcat(extractθ(gp, chains, ix)...))
end

function cond_latents(gp::LaplaceGP, θl, θc, x, y, z, x2)
	# Conditional latent distribution at the target x2

	# Evaluate the covariance function at θ
	xCx = zeros(size(x, 1), size(x, 1))
	covariance_function!(xCx, gp, θc, x)

	# Factorize
	L = cholesky(xCx, Val(false)).L
	@info "L" L

	# Laplace approximation for the latent posterior for the training points
	fwhat, ∇2lπ_fw = laplace_approx(gp, L, y, z, θl)
	@info "fhat" L*fwhat

	# Evaluate covariance between training and test points
	xCy = zeros(size(x2, 1), size(x, 1))
	covariance_function!(xCy, gp, θc, x2, x)
	@debug "xCy" xCy
	yCy = zeros(size(x2, 1), size(x2, 1))
	covariance_function!(yCy, gp, θc, x2)

	# Get the predicted distribution for the test point latents
	∇ll_f = zeros(length(y))
	∇2ll_f = zeros(length(y))
	fhat = L * fwhat
	ll = ∇2ll_f!(gp, ∇ll_f, ∇2ll_f, fhat, y, z, θl)

	@debug "∇2ll_f" ∇2ll_f

	# TODO: This was filled in with the naive formulas - this can be simplified
	# and the inv's can be removed
	ixCx = inv(xCx)
	μf2 = xCy * ixCx * fhat
	Σf2 = yCy .- xCy * (ixCx .- ixCx*inv(ixCx.+diagm(0=>∇2ll_f))*ixCx) * transpose(xCy)
	@debug "tr(Σf2)" diag(Σf2)

	return μf2, Σf2
end

function predict(gp::LaplaceGP, mcmc, x, y, z, x2, z2; quantiles=[], detail=15)
	# Get quadrature nodes and weights
	ghq_nodes, ghq_weights = gausshermite(detail)

	# Separate and transform parameters
	θl, θc = θ(gp, unrecord(gp, mcmc, size(mcmc.df,1)))

	# Get the latent distribution at the target points given the
	μf2, Σf2 = cond_latents(gp, θl, θc, x, y, z, x2)
	σ2f2 = diag(Σf2)

	# Gauss-Hermite quadrature to estimate predicted mean and variance at test points
	μ_pred = zeros(size(x2, 1))
	σ2_pred = zeros(size(x2, 1))
	Q_pred = zeros(size(x2, 1), length(quantiles))
	σf = sqrt.(σ2f)
	for i in 1:length(μf2)
		pred_nodes = [gp.datalik(μf[i] + nx * σf[i], z[i,:], θl) for nx in ghq_nodes]

		# Calculate mean
		μ_node = mean.(pred_nodes)
		μ = μ_pred[i] = sum(ghq_weights .* μ_node)

		# Calculate variance
		σ2_node = var.(pred_nodes)
		σ2 = σ2_pred[i] = sum(ghq_weights .* ((μ_node .- μ).^2 .+ σ2_node))

		# Quantiles
		σ = sqrt(σ2)
		for qi in eachindex(quantiles)
			# Initial guess
			minx, maxx = μ - σ, μ + σ
			target_q = quantiles[qi]
			q(x) = sum(ghq_weights .* cdf.(pred_nodes, x))

			# Expand search region until it encompases the target quantile
			while q(minx) > target_q
				minx, maxx = minx - 2 * (maxx - minx), minx
			end
			while q(maxx) < target_q
				minx, maxx = maxx, maxx + 2 * (maxx - minx)
			end

			# Bisecting search for the quantile
			for j in 1:36
				midx = (maxx + minx) / 2
				qmid = q(midx)
				if qmid < target_q
					minx = midx
				else
					maxx = midx
				end
			end

			# Round to integer if we're close enough
			qpred = (maxx + minx) / 2
			if qpred - round(qpred) < 1e-9
				qpred = round(qpred)
			end

			# And done..
			Q_pred[i,qi] = qpred
		end
	end

	return μ_pred, σ2_pred, Q_pred, μf2, σ2f2
end

function samplegp(gp::LaplaceGP, ϕ, x, y, z, x2, z2)
	# Separate and transform parameters
	θl, θc = θ(gp, ϕ)

	# Get the latent distribution at the target points given the
	μf2, Σf2 = cond_latents(gp, θl, θc, x, y, z, x2)

	# Sample the latents at the target points
	# Uses eigendecomposition to allow positive semi-definite covariance
	# matrices (could be optimized)
	E = eigen(Σf2)
	L = real.(E.vectors) * diagm(0 => sqrt.(max.(real.(E.values), 0.0)))
	f_samp = L * randn(length(μf2)) .+ μf2

	# Sample the data likelihood
	y_samp = rand.([gp.datalik(f_samp[i], z2[i,:], θl) for i in eachindex(f_samp)])

	return y_samp, f_samp
end

function logposterior(gp::LaplaceGP, ϕ, x, y, z)
	# Separate and transform parameters
	θl, θc, dθl_dϕl, dθc_dϕc = ∇θ(gp, ϕ)

	# Evaluate the covariance function at θ
	C = zeros(size(x, 1), size(x, 1))
	covariance_function!(C, gp, θc, x)

	# Factorize
	L = cholesky(C, Val(false)).L

	# Laplace approximation of the latent posterior
	fw, ∇2ll_fw, ll = laplace_approx(gp, L, y, z, θl)

	# Get the determinant of the covariance matrix from the cholesky
	# transform of the negative Hessian
	for i in 1:length(fw)
		for j in i:length(fw)
			# Negate and ensure symmetry
			∇2ll_fw[i,j] = ∇2ll_fw[j,i] = -(∇2ll_fw[i,j] + ∇2ll_fw[j,i]) / 2
		end
	end
	HL = cholesky!(∇2ll_fw, Val(false)).L
	ldetΣ = -sum(log.(diag(HL)))*2

	# Marginalize over the Laplace approximation of the latent posterior
	ll += (log(2*π) + ldetΣ) / 2

	# Log prior
	lp = sum(logpdf.(gp.θl_prior, θl) .+ log.(dθl_dϕl)) +
	     sum(logpdf.(gp.θc_prior, θc) .+ log.(dθc_dϕc))

	# Return log prior and log lik separately
	return lp, ll
end
