using CSV
using DataFrames
using DataFrames: transform, transform!
using Turing
using Distributions
using AbstractGPs, KernelFunctions
using LinearAlgebra
# using ReverseDiff
# using Zygote
# Turing.setadbackend(:forwarddiff)

ip1 = CSV.File("test/testin/input_pair_3206.tsv") |> DataFrame
ip1 = ip1[completecases(ip1), :]
pidmap = Dict(p=>i for (i,p) in enumerate(unique(ip1.PersonID)))
ip1.pid = [pidmap[p] for p in ip1.PersonID]

sekernel(alpha, rho) = 
  alpha^2 * KernelFunctions.transform(SEKernel(), sqrt(0.5)/rho)

function subjectcorrmat(subjects, t=T) where T
    n = length(subjects)
    scov = zeros(t, n, n)
    # if same person, 1.0, otherwise, 0.0
    for i in 1:n
        for j in 1:n
            subjects[i] == subjects[j] && (scov[i,j] = one(t))
        end
    end
    return scov
end

@model function GPmodel1(bug, diet, subj, tp, jitter=1e-6, T=Float64)
    nobs = length(bug)
    # assume zero mean
    mu = zeros(T, nobs)

    # subj-specific cov
    scov = subjectcorrmat(subj, T)

    # diet covariance
    dietcov = diet * diet'
    
    # variance contributed by diet slope
    b ~ Normal()
    dietcov *= b

    # gp priors
    sig2 ~ LogNormal(0, 1)
    alpha ~ LogNormal(0, 0.1)
    rho ~ LogNormal(0, 1)
    
    # gp specs
    kernel = sekernel(alpha, rho)  # covariance function
    K = kernelpdmat(kernel, reshape(tp, length(tp), 1), obsdim=1)  # cov matrix for time
    
    # remove cross-person covariance
    K = K .* scov
    # add outter product of diet * "slope" (cov caused by slope)
    K = K .+ dietcov
    
    # jitter for numerical stability
    K += LinearAlgebra.I * (sig2 + jitter)
    
    bug ~ MvNormal(mu, K)
end

gpm1 = GPmodel1(ip1.nutrient[1:150], ip1.bug[1:150], ip1.pid[1:150], ip1.Date[1:150])
@time r1 = sample(gpm1, HMC(0.1,20), 100)


@model function GPmodel2(bug, diet, subj, tp, jitter=1e-6)

    # diet coefficient prior
    b ~ Normal()
    # gp priors
    sig2 ~ LogNormal(0, 1)
    alpha ~ LogNormal(0, 0.1)
    rho ~ LogNormal(0, 1)


    for s in subj
        i = findall(==(s), subj)
        a ~ Normal()
        # lm specs
        mu = a .+ b * diet[i]
        
        # gp specs
        kernel = sekernel(alpha, rho)  # covariance function
        K = kernelpdmat(kernel, reshape(tp[i], length(tp[i]), 1), obsdim=1)  # cov matrix
        K += LinearAlgebra.I * (sig2 + jitter)
        
        bug[i] ~ Normal(mu, K)
    end
    
end
gpm2 = GPmodel2(ip1.nutrient[1:50], ip1.bug[1:50], ip1.pid[1:50], ip1.Date[1:50])
@time r2 = sample(gpm2, HMC(0.1,20), 100);


@model function mixedmodel1(bug, diet, subj)
    num_subjects = length(unique(subj))  # number of subjects
    
    intercept ~ Normal(0, 1)
    diet_effect ~ Normal(0, 1)
    subject_effects ~ filldist(Normal(0, 1), num_subjects)
    noise_sd ~ LogNormal(0, 1)  # eps_{s,j} ~ Normal(0, noise_sd)  

    m = intercept .+ diet * diet_effect + subject_effects[subj]
    bug .~ Normal.(m, noise_sd)  # you may need a different likelihood depending on the range of "bug".
end

mm1 = mixedmodel1(ip1.nutrient, ip1.bug, ip1.pid)
@time rmm1 = sample(mm1, HMC(0.1,20), 100)

@model function mixedmodel2(bug, subj)
    num_subjects = length(unique(subj))  # number of subjects
    
    intercept ~ Normal(0, 1)
    subject_effects ~ filldist(Normal(0, 1), num_subjects)
    noise_sd ~ LogNormal(0, 1)  # eps_{s,j} ~ Normal(0, noise_sd)  

    m = intercept .+ subject_effects[subj]
    bug .~ Normal.(m, noise_sd)  # you may need a different likelihood depending on the range of "bug".
end

mm2 = mixedmodel2(ip1.bug, ip1.pid)
@time rmm2 = sample(mm2, HMC(0.1,20), 100)



# log bayes
log2(mean(rmm1[:lp]) / mean(rmm2[:lp]))

rmm1[:diet_effect]

plot(rmm1)

############
# Plotting #
############

using StatsPlots
using Pipe
@pipe ip1 |> groupby(_, :pid) |> transform!(_, :pid => length => :nsamples)
filter!(:nsamples => ==(4), ip1)

@pipe ip1 |> groupby(_, :pid) |> 
        transform!(_, AsTable(:) => 
            (t-> (t.Date .- minimum(t.Date)) .+ 1) => :normDate)
ip1.normDate


@df ip1 plot(:normDate, :bug, group=:pid, legend=false,
             xlabel="Time", ylabel="Bug")

@df ip1 plot(:normDate, :nutrient, group=:pid, legend=false,
             xlabel="Time", ylabel="Diet")

@df ip1 scatter(:nutrient, :bug, legend=false,
                ylabel="Bug", xlabel="Diet")