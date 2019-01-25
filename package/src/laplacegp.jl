

struct LaplaceGP <: AbstractGP
    # The covariance function
    cf::Function
    # Covariance function hyperparameter transformation
    θc_link::Function
    # Prior on covariance function hyperparameters
    ϕc_prior::Vector{UnivariateDistribution}
    # Jitter
    jitter::Float64

    # Data Likelihood
    logdatalik::Function
    # Data likelihood hyperparameter transformation
    θl_link::Function
    # Prior on data likelihood hyperparameters
    ϕl_prior::Vector{UnivariateDistribution}
end

function covariance_function!(C::Matrix, gp::LaplaceGP, ϕc::Vector, x::Matrix)
end

function record(gp::LaplaceGP, chains::Chains, θ::Vector)
end

function predict(gp::LaplaceGP, θ::Vector, x::Matrix, x2::Matrix)
end

function ϕ(gp::LaplaceGP, θ::Vector)
    return gp.θl_link(θ[eachindex(gp.ϕl_prior)]),
        gp.θc_link(θ[length(gp.ϕl_prior) + eachindex(gp.ϕc_prior)])
end

function ∇2ll_lik_f!(∇ll, ∇2ll, logdatalik, f, ϕl)
	ll = 0.0
	for i = 1:length(f)
		@inbounds v = logdatalik(i, ForwardDiff.Dual(ForwardDiff.Dual(f[i], 1.0), ForwardDiff.Dual(1.0, 0.0)), ϕl)
		ll += ForwardDiff.value(ForwardDiff.value(v))
		@inbounds ∇ll[i] = ForwardDiff.partials(ForwardDiff.value(v), 1)
		@inbounds ∇2ll[i] = ForwardDiff.partials(ForwardDiff.partials(v, 1), 1)
	end

	return ll
end

function logposterior(gp::LaplaceGP, θ::Vector, x::Matrix, y::Vector, z::Vector)
    # Separate and transform parameters
    ϕl, ϕc = ϕ(gp, θ)

    # Evaluate the covariance function at θ
    C = zeros(size(x, 1))
    covariance_function!(C, gp, ϕc, x)

    # Factorize
    L = cholesky(C, Val(false)).L

	# Saddle-free Newton's method to find posterior mode for whitened latents fw
	fw = zeros(length(y))
    ∇2ll_fw = zeros(length(y), length(y))
	∇2ll_fw_temp = zeros(length(y), length(y))
	α = 1.0
	lltol = 0.01
	ll = -Inf
	maxN = 10

	for i in 1:maxN
		# Unwhiten
		mul!(f, L, fw)

		# Evaluate first and second derivatives of the likelihood at f
		ll2 = ∇2ll_lik_f!(∇ll_f, ∇2ll_f, gp.logdatalik, f, ϕl)
		@info "∇2ll_f", ∇2ll_f

		# Re-whiten Jacobian and Hessian
		mul!(∇ll_fw, transpose(L), ∇ll_f)
		# Slow formula: ∇2ll_fw = transpose(L) * Diagonal(∇2ll_f) * L
		mul!(∇2ll_fw_temp, diagm(0 => ∇2ll_f), L)
		mul!(∇2ll_fw, transpose(L), ∇2ll_fw_temp)

		# White noise prior on whitened latents
		ll2 -= sum(fw.^2) / 2.0
		∇ll_fw .-= fw
		for i in eachindex(fw)
			∇2ll_fw[i,i] -= 1.0
		end

		# Slow if we're actually getting worse
		α = ll2 < ll ? 0.5 : 1.0
		ll = ll2
		@info (i, ll)

		E = eigen(∇2ll_fw)
		@info "Eigs", E.values
		if any(E.values .> -1.0)
			# We're in an indefinite part of the parameter space, or where the
			# Hessian has dangerously low curvature
			# Apply a softabs transform on eigenvalues to avoid saddle
			# points/minima and to regularize small values in the Hessian
			λ = E.values .* coth.(E.values) # softabs with minimum value 1
			∇2ll_fw_reg .= E.vectors * diagm(0 => -λ) * transpose(E.vectors)
			@warn "Applying softabs" λ
		else
			# We're in a safely-positive definite part of the parameter space
			∇2ll_fw_reg .= ∇2ll_fw

			# How much of a gain do we expect?
			#ldiv!(dfw, ∇2ll_fw, ∇ll_fw)
			dfw .= ∇2ll_fw \ ∇ll_fw
			llgain = -transpose(dfw) * ∇2ll_fw * dfw
			@info "Gain" llgain
			llgain < lltol && break # stop when we gain too little from another step
		end

		# Newton step
		if i < maxN
			fw .-= α * (∇2ll_fw_reg \ ∇ll_fw)
		else
			@warn "maxN reached in Newton steps"
		end
	end

	# Get the determinant of the covariance matrix from the cholesky
	# transform of the negative Hessian
	for i in 1:Nf
		for j in i:Nf
			# Negate and ensure symmetry
			∇2ll_fw[i,j] = ∇2ll_fw[j,i] = -(∇2ll_fw[i,j] + ∇2ll_fw[j,i]) / 2
		end
	end
	HL = cholesky!(∇2ll_fw, Val(false)).L
	ldetΣ = -sum(log.(diag(HL)))*2

	# Marginalize over the Laplace approximation of the latent posterior
	ll += (log(2*π) + ldetΣ) / 2

    # Log prior
    lp = sum(logpdf.(ϕl_prior, ϕl)) + sum(logpdf.(ϕc_prior, ϕc))

    # Log posterior
    return lp + ll
end
