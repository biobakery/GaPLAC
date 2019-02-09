

struct LaplaceGP <: AbstractGP
    # The covariance function
	# Signature: (x_i::Vector, x_j::Vector, self::Boolean, θc::Vector) -> Real
    cf::Function
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

function θ(gp::LaplaceGP, ϕ::Vector)
    return gp.θl_link(ϕ[eachindex(gp.θl_prior)]),
        gp.θc_link(ϕ[length(gp.θl_prior) + eachindex(gp.θc_prior)])
end

function ∇θ(gp::LaplaceGP, ϕ::Vector)
	θld, θcd = θ([ForwardDiff.Dual(ϕi, 1.) for ϕi in ϕ])
    return value.(θld), value.(θcd), partials.(θld,1), partials.(θcd,1)
end


function covariance_function!(C::Matrix, gp::LaplaceGP, θc::Vector, x::Matrix)
	for i in 1:size(x,1)
		for j in i:size(x,1)
			C[i,j] = C[j,i] = gp.cf(x[i,:], x[j,:], i==j)
		end
	end
end

function covariance_function!(C::Matrix, gp::LaplaceGP, θc::Vector, x1::Matrix, x2::Matrix)
	for i in 1:size(x1,1)
		for j in 1:size(x2,1)
			C[i,j] = gp.cf(x1[i,:], x2[j,:], false)
		end
	end
end

function ∇2ll_f!(gp::LaplaceGP, ∇ll::Vector, ∇2ll::Vector,
		f::Vector, z::Vector, y::Vector, θl::Vector)
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

function laplace_approx(gp::LaplaceGP, L::LowerTriangular, y::Vector, θl::Vector;
	lπtol = 0.01, maxN = 10)

	# Saddle-free Newton's method to find posterior mode for whitened latents fw
	fw = zeros(length(y))
	∇2ll_fw = zeros(length(y), length(y))
	∇2ll_fw_temp = zeros(length(y), length(y))
	α = 1.0
	lπ = -Inf

	for i in 1:maxN
		# Unwhiten
		mul!(f, L, fw)

		# Evaluate first and second derivatives of the likelihood at f
		ll = ∇2ll_lik_f!(∇ll_f, ∇2ll_f, gp.logdatalik, f, θl)
		@info "∇2ll_f", ∇2ll_f

		# Re-whiten Jacobian and Hessian
		mul!(∇ll_fw, transpose(L), ∇ll_f)
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
		@info (i, lπ)

		E = eigen(∇2lπ_fw)
		@info "Eigs", E.values
		if any(E.values .> -1.0)
			# We're in an indefinite part of the parameter space, or where the
			# Hessian has dangerously low curvature
			# Apply a softabs transform on eigenvalues to avoid saddle
			# points/minima and to regularize small values in the Hessian
			λ = E.values .* coth.(E.values) # softabs with minimum value 1
			∇2lπ_fw_reg .= E.vectors * diagm(0 => -λ) * transpose(E.vectors)
			@warn "Applying softabs" λ
		else
			# We're in a safely-positive definite part of the parameter space
			∇2lπ_fw_reg .= ∇2lπ_fw

			# How much of a gain do we expect?
			#ldiv!(dfw, ∇2ll_fw, ∇ll_fw)
			dfw .= ∇2lπ_fw \ ∇lπ_fw
			lπgain = -transpose(dfw) * ∇2lπ_fw * dfw
			@info "Gain" lπgain
			lπgain < lπtol && break # stop when we gain too little from another step
		end

		# Newton step
		if i < maxN
			fw .-= α * (∇2ll_fw_reg \ ∇ll_fw)
		else
			@warn "maxN reached in Newton steps"
		end
	end

	return fw, ∇2lπ_fw, lπ
end

function record(gp::LaplaceGP, chains::Chains, ϕ::Vector)
	# Separate and transform parameters
    θl, θc = θ(gp, ϕ)
end

function predict(gp::LaplaceGP, ϕ::Vector, x::Matrix, y::Vector, x2::Matrix;
		quantiles::Vector)
	# Separate and transform parameters
    θl, θc = θ(gp, ϕ)

	# Evaluate the covariance function at θ
    xCx = zeros(size(x, 1), size(x, 1))
    covariance_function!(xCx, gp, θc, x)

	# Factorize
	L = cholesky(xCx, )

	# Laplace approximation for the latent posterior for the training points
	fw, ∇2lπ_fw = laplace_approx(gp, L, y, θl)

	# Evaluate covariance between training and test points
	xCy = zeros(size(x, 1), size(x2, 1))
    covariance_function!(xCy, gp, θc, x, x2)
	yCy = zeros(size(x2, 1), size(x2, 1))
    covariance_function!(yCy, gp, θc, x2, x2)

	# Get the predicted distribution for the test point latents
	???
	μf = ???
	σ2f = diag(???)

	# Gaussian quadrature to estimate predicted mean and variance at test points
	μ_pred = zeros(size(x2, 1))
	σ2_pred = zeros(size(x2, 1))
	σf = sqrt.(σ2f)
	for i in 1:size(x2, 1)
		for <quadrature loop>
			# Make the prediction
			pred = gp.datalik(μf[i] + ??? * σf[i], z[i,:], θl)

			# Mix the prediction
			mean(pred) # TODO
			var(pred)
		end

		# TODO: Find the quantiles
	end

	return μ_pred, σ2_pred, Q_pred, μf, σ2f
end

function logposterior(gp::LaplaceGP, ϕ::Vector, x::Matrix, y::Vector, z::Vector)
    # Separate and transform parameters
    θl, θc, dθl_dϕl, dθc_dϕc = ∇θ(gp, ϕ)

    # Evaluate the covariance function at θ
    C = zeros(size(x, 1), size(x, 1))
    covariance_function!(C, gp, θc, x)

    # Factorize
    L = cholesky(C, Val(false)).L

	# Laplace approximation of the latent posterior
	fw, ∇2ll_fw, ll = laplace_approx(gp, L, θl)

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
    lp = sum(logpdf.(θl_prior, θl) .+ log.(dθl_dϕl)) +
	     sum(logpdf.(θc_prior, θc) .+ log.(dθc_dϕc))

    # Log posterior
    return lp + ll
end
