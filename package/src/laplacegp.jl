
using LinearAlgebra

function θ(gp::LaplaceGP, ϕ)
    return gp.θl_link(ϕ[1:length(gp.θl_prior)]),
        gp.θc_link(ϕ[length(gp.θl_prior) .+ (1:length(gp.θc_prior))])
end

function ∇θ(gp::LaplaceGP, ϕ)
	θld, θcd = θ([ForwardDiff.Dual(ϕi, 1.) for ϕi in ϕ])
    return value.(θld), value.(θcd), partials.(θld,1), partials.(θcd,1)
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
	∇ll_f, ∇ll_fw = zeros(length(y)), zeros(length(y))
	∇lπ_fw = zeros(length(y))
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
		@info "∇2ll_f" ∇2ll_f

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

function record!(gp::LaplaceGP, chains::Chains, ϕ)
	# Separate and transform parameters
    θl, θc = θ(gp, ϕ)

	for i in eachindex(θl)
		record!(chains, Symbol(θl_names[i]), θl[i])
	end
	for i in eachindex(θc)
		record!(chains, Symbol(θc_names[i]), θc[i])
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
		maxv[θ_mid .> θ_target] .= midv[θ_mid .> θ_target]
		minv[θ_mid .< θ_target] .= midv[θ_mid .< θ_target]
		n += 1
	end

	return maxv
end

function unrecord(gp::LaplaceGP, chains::Chains, ix::Int)
	# Build the target θ vector
	θ_target = vcat([getrecords(chains, Symbol(name))[ix] for name in gp.θl_names],
		[getrecords(chains, Symbol(name))[ix] for name in gp.θc_names])

	return invlink(gp, θ_target)
end

function cond_latents(gp::LaplaceGP, θl, θc, x, y, z, x2)
	# Conditional latent distribution at the target x2

	# Evaluate the covariance function at θ
    xCx = zeros(size(x, 1), size(x, 1))
    covariance_function!(xCx, gp, θc, x)

	# Factorize
	L = cholesky(xCx, Val(false)).L

	# Laplace approximation for the latent posterior for the training points
	fwhat, ∇2lπ_fw = laplace_approx(gp, L, y, z, θl)

	# Evaluate covariance between training and test points
	xCy = zeros(size(x2, 1), size(x, 1))
    covariance_function!(xCy, gp, θc, x2, x)
	yCy = zeros(size(x2, 1), size(x2, 1))
    covariance_function!(yCy, gp, θc, x2)

	# Get the predicted distribution for the test point latents
	∇ll_f = zeros(length(y))
	∇2ll_f = zeros(length(y))
	ll = ∇2ll_f!(gp, ∇ll_f, ∇2ll_f, L*fwhat, y, z, θl)

	μf2 = xCy * ∇ll_f
	Σf2 = yCy .- xCy * (∇2lπ_fw \ transpose(xCy))

	return μf2, Σf2
end

function predict(gp::LaplaceGP, ϕ, x, y, z, x2; quantiles)

	# Separate and transform parameters
    θl, θc = θ(gp, ϕ)

	# Get the latent distribution at the target points given the
	μf2, Σf2 = cond_latents(gp, θl, θc, x, y, z, x2)

	# Gaussian quadrature to estimate predicted mean and variance at test points
	μ_pred = zeros(size(x2, 1))
	σ2_pred = zeros(size(x2, 1))
	σf = sqrt.(σ2f)
	for i in 1:size(x2, 1)
		begin#for <quadrature loop>
			# Make the prediction
			#pred = gp.datalik(μf[i] + ??? * σf[i], z[i,:], θl)

			# Mix the prediction
			mean(pred) # TODO
			var(pred)
		end

		# TODO: Find the quantiles
	end

	return μ_pred, σ2_pred, Q_pred, μf, σ2f
end

function samplegp(gp::LaplaceGP, ϕ, x, y, z, x2, z2)
	# Separate and transform parameters
    θl, θc = θ(gp, ϕ)

	# Get the latent distribution at the target points given the
	μf2, Σf2 = cond_latents(gp, θl, θc, x, y, z, x2)

	# Sample the latents at the target points
	f_samp = rand(MvNormal(μf2, Σf2))

	# Sample the data likelihood
	y_samp = rand.([gp.datalik(f_samp[i], z2[i], θl) for i in eachindex(z2)])

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

    # Return log prior and log lik separately
    return lp, ll
end
