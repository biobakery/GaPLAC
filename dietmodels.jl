using CSV
using DataFrames
using DataFrames: transform, transform!
using Turing
using Distributions
using AbstractGPs, KernelFunctions
using LinearAlgebra
using Chain
using Random
using Logging, LoggingExtras

# global_logger(SimpleLogger(stderr, Logging.Debug))
# using ReverseDiff
# using Zygote
# Turing.setadbackend(:forwarddiff)

function subjectcorrmat(subjects, t=T) where T
    n = length(subjects)
    scov = zeros(t, n, n)
    # if same person, true, otherwise, false
    for i in 1:n, j in 1:n
        subjects[i] == subjects[j] && (scov[i,j] = one(t))
    end
    return scov
end

sekernel(α, ρ) = α^2 * KernelFunctions.transform(SEKernel(), 1/ρ)
maternkernel(α, ρ) = α^2 * KernelFunctions.transform(Matern52Kernel(), 1/ρ)
linearkernel(m, c) = KernelFunctions.transform(LinearKernel(c=c), m)

function get_cov(tp, σ, α, ρ, κ, jitter)
    @debug size(tp)
    n_obs = length(tp)
    # subject-specific covariance
    # subject_cov = subjectcorrmat(subject, Float64)
    # cov matrix for time, size n_obs x n_obs, obsdim=1 specifies that the matrix (tp?) is in the form of number of samples x number of features
    # Distance matrix.
    K = kernelpdmat(κ(α, ρ), reshape(tp, length(tp), 1), obsdim=1)
    # K = kernelpdmat(Matern52Kernel(), reshape(tp, length(tp), 1), obsdim=1)
    # for some reason if we don't remove the subject_cov, we get Float64 errors, probably having to do with the weird format of K after kernelpdmat. This converts it to a regular Float64 Matrix.
    K = K * Matrix{Float64}(LinearAlgebra.I, n_obs, n_obs)
    # remove cross-subject covariance
    # K = K .* subject_cov
     # jitter for numerical stability
    return K + LinearAlgebra.I * (σ + jitter)
end

function log2bayes(s1, s2)
    l2b = log2(mean(map(x-> BigFloat(2) ^ x, s1[:lp])) / mean(map(x-> BigFloat(2) ^ x, s2[:lp])))
    m1 = mean(s1[:lp])
    m2 = mean(s2[:lp])
    return l2b, m1, m2
end

@model function GPmodel1(bug, tp, subject_ids, kernel=sekernel, jitter=1e-6, T=Float64)
    # Priors.
    alpha ~ filldist(LogNormal(0.0, 0.1), length(subject_ids))
    rho ~ filldist(LogNormal(0.0, 1.0), length(subject_ids))
    sigma ~ LogNormal(0.0, 1.0)
    
    for (index, subject_id) in enumerate(subject_ids)

        n_obs = length(bug[index])
        mu = zeros(Float64, n_obs)

        K = get_cov(tp[index], sigma, alpha[index], rho[index], kernel, jitter)
        # println(size(K))

        bug[index] ~ MvNormal(mu, K)
    end
end

@model function GPmodel2(bug, diet, tp, subject_ids, kernel=sekernel, jitter=1e-6, T=Float64)
    # Priors.
    alpha ~ filldist(LogNormal(0.0, 0.1), length(subject_ids))
    rho ~ filldist(LogNormal(0.0, 1.0), length(subject_ids))
    sigma ~ LogNormal(0.0, 1.0)

    # diet
    m ~ LogNormal(0, 1)
    c ~ LogNormal(0, 1)
    
    
    for (index, subject_id) in enumerate(subject_ids)

        n_obs = length(bug[index])
        mu = zeros(Float64, n_obs)

        K_time = get_cov(tp[index], sigma, alpha[index], rho[index], kernel, jitter)
        K_diet = get_cov(diet[index], sigma, m, c, linearkernel, jitter)
        # println(size(K))

        bug[index] ~ MvNormal(mu, K_time .+ K_diet)
    end
end
