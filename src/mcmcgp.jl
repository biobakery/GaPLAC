
function mcmcgp(gp::AbstractGP, x, y, z, θ_start, Nsamples)

	# Setup chains
	chain = Chains()

	# Actual chains
	startTime = time_ns()
	mcmc = UnivariateSliceMCMC(ϕ -> begin
		lp, ll = logposterior(gp, ϕ, x, y, z)
		lp + ll
	end, repeat([0.5], inner=length(θ_start)), 1)
	stats = MCMCStats()
	ϕ = invlink(gp, θ_start)
	ϕsamp = copy(ϕ)
	lπ = mcmc.f(ϕ)

	# TODO: Implement annealing

	for sampi = 0:Nsamples
		if sampi >= 1
			try
				# Generate one sample
				copyto!(ϕsamp, ϕ)
				lπ = mcmcsample!(ϕsamp, mcmc, stats)
				copyto!(ϕ, ϕsamp)
			catch err
				@warn(@sprintf("MCMC sample generation failed: %s", string(err)))
				continue
			end
		end

		# Record sample
		newsample!(chain)
		record!(gp, chain, ϕ)
		record!(chain, :lπ, lπ)

		# Progress
		if sampi == 1 || sampi == 10 || sampi > 0 && floor(10 * (sampi - 1) / Nsamples) < floor(10 * sampi / Nsamples)
			elapsed = (time_ns() - startTime) / 1e9
			if sampi == 1
				remaining = (Nsamples - sampi) * ((elapsed - stats.total_time_s)/2 + stats.total_time_s)
			else
				remaining = (Nsamples - sampi) * (elapsed / sampi)
			end

			@info "MCMC: $(sampi)/$Nsamples samples ($(@sprintf("%.0f", elapsed))s; ~$(@sprintf("%.0f", remaining))s left) lπ=$(@sprintf("%.2f", lπ))"
		end
	end

	elapsed = (time_ns() - startTime) / 1e9
	sampsPerSec = Nsamples / elapsed
	@info "$(@sprintf("%.2f", sampsPerSec)) samples/s; time: $(@sprintf("%.1f", 100*stats.total_time_s/elapsed))% samp, $(@sprintf("%.1f", 100*(elapsed-stats.total_time_s)/elapsed))% other"
	@info "Acceptance rate: $(@sprintf("%d / %d (%.1f%%)", stats.accepts, stats.accepts + stats.rejects, 100*stats.accepts/(stats.accepts + stats.rejects)))"

	return chain
end
