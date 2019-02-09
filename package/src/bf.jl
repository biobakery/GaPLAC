
function hmlπ(c::Chains)
    # Estimate the harmonic mean of the unnormalized log posterior
    # Note: log base is kept as the natural log
    lπ = getrecords(c, :lπ)

    # Shift the range so we can exponentiate
    a = maximum(lπ)
    reg_lπ = max(log(.5 / length(lπ)), lπ .- a) # truncate unreasonably rare events

    # Regularized harmonic mean of lπ
    reg_hmlπ = -log(mean(exp.(.-reg_lπ)))

    # Unshift the scale
    lπ + a
end

function log2_bayes_factor(c1::Chains, c2::Chains)
    # Log2 Bayes factor between MCMC chains drawn from two models
    (hmlπ(c1) - hmlπ(c2)) / log(2)
end
