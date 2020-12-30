using CSV
using DataFrames
using DataFrames: transform, transform!
using Turing
using Distributions
using AbstractGPs, KernelFunctions
using LinearAlgebra
using Chain
using Random

# using ReverseDiff
# using Zygote
# Turing.setadbackend(:forwarddiff)

function subjectcorrmat(subjects)
    n = length(subjects)
    scov = zeros(Bool, n, n)
    # if same person, true, otherwise, false
    for i in 1:n, j in 1:n
        subjects[i] == subjects[j] && (scov[i,j] = true)
    end
    return scov
end

@model function GPmodel1(bug, diet, subjcorr, tp, jitter=1e-6, T=Float64)
    nobs = length(bug)
    # assume zero mean
    μ = zeros(T, nobs)
    
    # diet covariance
    c ~ LogNormal(0, 1)
    K_diet = kernelmatrix(LinearKernel(c=c), reshape(diet, length(diet),1), obsdim=1)
    
    # time covariance
    K_time = kernelmatrix(Matern52Kernel(), reshape(tp, length(tp), 1), obsdim=1)  # cov matrix for time
    # remove cross-person covariance
    K_time = K_time .* subjcorr
    
    # jitter for numerical stability
    σ2 ~ LogNormal(0, 1)
    # K_diet += LinearAlgebra.I * (σ2 + jitter)
    K_time += LinearAlgebra.I * (σ2 + jitter)

    bug ~ MvNormal(μ, K_time .+ K_diet)
end

@model function GPmodel2(bug, subjcorr, tp, jitter=1e-6, T=Float64)
    nobs = length(bug)
    # assume zero mean
    μ = zeros(T, nobs)
        
    # time covariance
    K_time = kernelmatrix(Matern52Kernel(), reshape(tp, length(tp), 1), obsdim=1)  # cov matrix for time
    # remove cross-person covariance
    K_time = K_time .* subjcorr
    
    # jitter for numerical stability
    σ2 ~ LogNormal(0, 1)
    # K_diet += LinearAlgebra.I * (σ2 + jitter)
    K_time += LinearAlgebra.I * (σ2 + jitter)

    bug ~ MvNormal(μ, K_time)
end

log2bayes(s1, s2) = log2(mean(map(x-> 2. ^ x, s1[:lp])) / mean(map(x-> 2. ^ x, s2[:lp])))