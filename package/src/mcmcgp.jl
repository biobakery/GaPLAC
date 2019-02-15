
function mcmcgp(gp::AbstractGP, x, y, z, ϕ, Nsamples)

	# Setup chains
	chain = Chains()

	# Actual chains
	startTime = time_ns()
	# TODO: Implement annealing
	mcmc = UnivariateSliceMCMC(ϕ -> begin
		lp, ll = logposterior(gp, ϕ, x, y, z)
		lp + ll
	end, 0.5, 1)
	stats = MCMCStats()
	logpost = 0
	ϕsamp = copy(ϕ)

	for sampi = 0:Nsamples
		if sampi >= 1
			try
				# Generate one sample
				copyto!(ϕsamp, ϕ)
				lπ = sample!(ϕsamp, mcmc, stats)
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
				remaining = (Nsamples - sampi) * ((elapsed-sampTime/1e9)/2 + sampTime/1e9)
			else
				remaining = (Nsamples - sampi) * (elapsed / sampi)
			end

			@info "MCMC: $(sampi)/$Nsamples samples ($(@sprintf("%.0f", elapsed))s; ~$(@sprintf("%.0f", remaining))s left) lpost=$(@sprintf("%.2f", logpost))"
		end
	end

	elapsed = (time_ns() - startTime) / 1e9
	sampsPerSec = Nsamples / elapsed
	sampTime /= 1e9
	@info "$(@sprintf("%.2f", sampsPerSec)) samples/s; time: $(@sprintf("%.1f", 100*sampTime/elapsed))% samp, $(@sprintf("%.1f", 100*(elapsed-sampTime)/elapsed))% other"
	@info "Acceptance rate: $(@sprintf("%d / %d (%.1f%%)", mcmc.accepts, mcmc.accepts + mcmc.rejects, 100*mcmc.accepts/(mcmc.accepts + mcmc.rejects)))"

	return chain
end
