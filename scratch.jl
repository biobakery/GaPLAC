using CSV
using DataFrames
using DataFrames: transform, transform!
using Turing
using Distributions
using AbstractGPs, KernelFunctions
using LinearAlgebra
using StatsPlots
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

## -- Fake Data -- ##

function run_sim(nsub, dietdist=Normal())
    subj = repeat(1:nsub, inner=4)
    tp = repeat([1,2,9,10], outer=nsub) 
    bugind = let base = randn(4*nsub)
        for (i, b) in enumerate(base)
            if isodd(i % 4)
                base[i] = b + base[i+1]
            end
        end
        base
    end
    diet = rand(dietdist, 4*nsub)
    bugdiet = bugind .+ 1 .* diet # diet effect
    bugdiet2 = bugind .+ 2 .* diet # bigger diet effect
    sc = subjectcorrmat(subj)

    ind_gp1 = GPmodel1(bugind, diet, sc, tp)
    ind_samp1 = sample(ind_gp1, HMC(0.1,20), 100)

    ind_gp2 = GPmodel2(bugind, sc, tp)
    ind_samp2 = sample(ind_gp2, HMC(0.1,20), 100)

    ind_l2b = log2bayes(ind_samp1, ind_samp2)

    diet_gp1 = GPmodel1(bugdiet, diet, sc, tp)
    diet_samp1 = sample(diet_gp1, HMC(0.1,20), 100)

    diet_gp2 = GPmodel2(bugdiet, sc, tp)
    diet_samp2 = sample(diet_gp2, HMC(0.1,20), 100)

    diet_l2b = log2bayes(diet_samp1, diet_samp2)

    diet2_gp1 = GPmodel1(bugdiet2, diet, sc, tp)
    diet2_samp1 = sample(diet2_gp1, HMC(0.1,20), 100)

    diet2_gp2 = GPmodel2(bugdiet2, sc, tp)
    diet2_samp2 = sample(diet2_gp2, HMC(0.1,20), 100)

    diet2_l2b = log2bayes(diet2_samp1, diet2_samp2)
    
    return (ind=ind_l2b, diet1=diet_l2b, diet2=diet2_l2b)
end
 
run_sim(4)
df = DataFrame([run_sim(4) for _ in 100])

## -- Real Data -- ##

ip1 = CSV.File("test/testin/input_pair_1609.tsv") |> DataFrame
ip1 = ip1[completecases(ip1), :]
pidmap = Dict(p=>i for (i,p) in enumerate(unique(ip1.PersonID)))
ip1.pid = [pidmap[p] for p in ip1.PersonID]
ip1.datemod = [d+rand(Normal(0, 0.2)) for d in ip1.Date]
ip1.dietmod = [d+rand(Normal(0, 0.2)) for d in ip1.nutrient]
ip1.bugmod = [b+rand(Normal(0, 0.2)) for b in ip1.bug]
@chain ip1 begin
    groupby(:pid)
    transform!(:pid => length => :nsamples)
end

scorr = subjectcorrmat(ip1.pid)
gpm1 = GPmodel1(ip1.bugmod, ip1.nutrient, scorr, ip1.datemod)
@time r1 = sample(gpm1, HMC(0.1,20), 100)

plot(r1)

gpm2 = GPmodel2(ip1.bug, scorr, ip1.datemod)
@time r2 = sample(gpm2, HMC(0.1,20), 100);

# log bayes
r1[:lp]
log2(mean(map(x-> 2. ^ x, r1[:lp])) / mean(map(x-> 2. ^ x, r2[:lp])))

log2(mean(r1[:lp] .^2) / mean(r2[:lp] .^2))


@model function mixedmodel1(bug, diet, subj)
    num_subjects = length(unique(subj))  # number of subjects
    
    intercept ~ Normal(0, 1)
    diet_effect ~ Normal(0, 1)
    subject_effects ~ filldist(Normal(0, 1), num_subjects)
    noise_sd ~ LogNormal(0, 1)  # eps_{s,j} ~ Normal(0, noise_sd)  

    m = intercept .+ diet * diet_effect + subject_effects[subj]
    bug .~ Normal.(m, noise_sd)
end

mm1 = mixedmodel1(ip1.nutrient, ip1.bug, ip1.pid)
@time rmm1 = sample(mm1, HMC(0.1,20), 100)

# log bayes GP vs mixed model
log2(mean(r1[:lp] .^2) / mean(rmm1[:lp] .^2))

@model function mixedmodel2(bug, subj)
    num_subjects = length(unique(subj))  # number of subjects
    
    intercept ~ Normal(0, 1)
    subject_effects ~ filldist(Normal(0, 1), num_subjects)
    noise_sd ~ LogNormal(0, 1)  # eps_{s,j} ~ Normal(0, noise_sd)  

    m = intercept .+ subject_effects[subj]
    bug .~ Normal.(m, noise_sd)
end

mm2 = mixedmodel2(ip1.bug, ip1.pid)
@time rmm2 = sample(mm2, HMC(0.1,20), 100)


# log bayes
log2(mean(rmm1[:lp] .^2) / mean(rmm2[:lp] .^2))


ip2 = CSV.File("test/testin/input_pair_109.tsv") |> DataFrame
ip2 = ip2[completecases(ip2), :]
pidmap = Dict(p=>i for (i,p) in enumerate(unique(ip2.PersonID)))
ip2.pid = [pidmap[p] for p in ip2.PersonID]
ip2.datemod = [d+rand(Normal(0, 0.1)) for d in ip2.Date]

scorr2 = subjectcorrmat(ip2.pid)
gpm12 = GPmodel1(ip2.bug, ip2.nutrient, scorr2, ip2.datemod)
@time r12 = sample(gpm12, HMC(0.1,20), 100)
gpm22 = GPmodel2(ip2.bug, ip2.pid, ip2.datemod)
@time r22 = sample(gpm22, HMC(0.1,20), 100);

# log bayes
log2(mean(r12[:lp] .^2) / mean(r22[:lp] .^2))


################
# Benchmarking #
################
using BenchmarkTools

@model function exteriorkernel(bug, K_time, K_diet, jitter=1e-6, T=Float64)
    nobs = length(bug)
    # assume zero mean
    μ = zeros(T, nobs)
      
    # jitter for numerical stability
    σ2 ~ LogNormal(0, 1)
    K_time += LinearAlgebra.I * (σ2 + jitter)

    bug ~ MvNormal(μ, K_time .+ K_diet)
end

@model function exteriortime(bug, diet, K_time, jitter=1e-6, T=Float64)
    nobs = length(bug)
    # assume zero mean
    μ = zeros(T, nobs)
    
    # diet covariance
    c ~ LogNormal(0, 1)
    K_diet = kernelmatrix(LinearKernel(c=c), reshape(diet, length(diet),1), obsdim=1)
    
    # jitter for numerical stability
    σ2 ~ LogNormal(0, 1)
    K_time += LinearAlgebra.I * (σ2 + jitter)

    bug ~ MvNormal(μ, K_time .+ K_diet)
end

s = ip1.pid[1:200]
b = ip1.bug[1:200]   
d = ip1.nutrient[1:200]    
t = ip1.datemod[1:200]    
sc = subjectcorrmat(s)
   
K_diet = kernelmatrix(LinearKernel(), reshape(d, length(d),1), obsdim=1)
K_time = kernelmatrix(Matern52Kernel(), reshape(t, length(t), 1), obsdim=1)  # cov matrix for time
# remove cross-person covariance
K_time = K_time .* sc

gp1 = GPmodel1(b, d, sc, t)
gp2 = exteriorkernel(b, K_time, K_diet)
gp3 = exteriortime(b, d, K_time)

@btime sample($gp1, HMC(0.1,20), 100)
@btime sample($gp2, HMC(0.1,20), 100)
@btime sample($gp3, HMC(0.1,20), 100)

############
# Plotting #
############

using StatsPlots

@chain ip1 begin
    groupby(:pid)
    transform!(:pid => length => :nsamples)
end

@chain ip1 begin
    groupby(:pid)
    transform!(
        AsTable(:) => 
        (t-> (t.Date .- minimum(t.Date)) .+ 1) => 
        :normDate)
end

@chain ip1 begin
    groupby(:pid)
    transform!(
        AsTable(:) => 
        (t-> (t.datemod .- minimum(t.datemod)) .+ 1) => 
        :normDatemod)
end

@df ip1 scatter(:normDate, :bug, color=:gray,
             xlabel="Time", ylabel="Bug", label="Date")
@df ip1 scatter!(:normDatemod, :bug, color=:blue, alpha=0.2,
            xlabel="Time", ylabel="Bug", label="Date + rand")

@df filter(:normDate => <(20), ip1) scatter(:normDate, :bug, color=:gray,
                                            xlabel="Time", ylabel="Bug", label="Date")
@df filter(:normDate => <(20), ip1) scatter!(:normDatemod, :bug, color=:blue, alpha=0.2,
                                            xlabel="Time", ylabel="Bug", label="Date + rand")


@df ip1 scatter(:bug, :nutrient, color=:gray,
            xlabel="Bug", ylabel="Nutrient", label="Diet")
@df ip1 scatter!(:bug, :dietmod, color=:blue, alpha=0.2,
            xlabel="Bug", ylabel="Nutrient", label="Diet + rand")
                               