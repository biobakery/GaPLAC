using KernelFunctions: Distances
using GaPLAC
using CairoMakie
using AbstractGPs
using KernelFunctions
using AbstractGPsMakie
using Distances

using CSV
using DataFrames

struct IdentityKernal <: KernelFunctions.SimpleKernel end
KernelFunctions.kappa(::IdentityKernal, d2::Real) = d2 > 0 ? 0 : 1
KernelFunctions.metric(::IdentityKernal) = Euclidean()

data = CSV.read("test/testin/input_pair_3206.tsv", DataFrame)
filter!(:Date=> !ismissing, data)
data.Date = disallowmissing(data.Date)
gdf = groupby(data, :PersonID)

# Give each subject a random offset from 0
data = transform(gdf, :Date=> (d-> d .- minimum(d) .+ rand()) => :DateOffset)

fig = Figure(resolution=(1200, 600))
ax1 = Axis(fig[1,1], xlabel="bug", ylabel="nutrient")
ax2 = Axis(fig[1,2], xlabel="date (offset)", ylabel="bug")
ax3 = Axis(fig[1,3], xlabel="date (offset)", ylabel="nutrient")

scatter!(ax1, data.bug, data.nutrient)
scatter!(ax2, data.DateOffset, data.bug)
scatter!(ax3, data.DateOffset, data.nutrient)

fig

k = Matern52Kernel()
gp = GP(k)

fx = gp(df.bug, 0.01)
logpdf(fx, df.nutrient)

p_fx = posterior(fx, df.nutrient)

plot!(ax1, -1:0.1:3, p_fx; bandscale=3, color=(:blue, 0.3))
fig


param_1, param_2, param_3, param_4, noise = ones(5)

k_nutrient = LinearKernel() ∘ ScaleTransform(param_1) # kernel for nutrient
# k_obese = transform(LinearKernel(), param_2) # kernel for obese
k_t = SqExponentialKernel() ∘ ScaleTransform(param_3) # kernel for time
k_sub = IdentityKernal() ∘ ScaleTransform(param_4) # kernel for subject

k = k_nutrient ⊗ k_t ⊗ k_sub # Collect all the kernels to make them act dimension wise
# Here we create a the prior based on the kernel and the data
priorgp = AbstractGPs.FiniteGP(GP(k), hcat(data.nutrient, data.DateOffset, data.PersonID) , noise, obsdim = 1) # Not sure about the dimensionality here
# We finally compute the posterior given the y observations
posteriorgp = posterior(priorgp, data.bug)

logpdf(priorgp, data.bug)

k2 = k_t ⊗ k_sub
priorgp2 = AbstractGPs.FiniteGP(GP(k2), hcat(data.DateOffset, data.PersonID) , noise, obsdim = 1) # Not sure about the dimensionality here
posteriorgp2 = posterior(priorgp2, data.bug)

logpdf(priorgp2, data.bug)

logpdf(posteriorgp(hcat(data.nutrient, data.DateOffset, data.PersonID), obsdim=1), data.bug)

using StatsBase


mc = mean_and_cov(posteriorgp, RowVecs(hcat(data.nutrient, data.DateOffset, data.PersonID)))

## Example for zulip

using AbstractGPs
using KernelFunctions
using CairoMakie
using AbstractGPsMakie

k1 = LinearKernel()
k2 = SqExponentialKernel()

gp = GP(k1 ⊗ k2)
x = randn(10)
t = 10 .* rand(10)

gpxt = AbstractGPs.FiniteGP(gp, RowVecs(hcat(x, t)), 0.1)

y = rand(gpxt)

fig = Figure(resolution=(1500, 600))
ax1 = Axis(fig[1,1], xlabel="time", ylabel="x")
ax2 = Axis(fig[1,2], xlabel="time", ylabel="y")
ax3 = Axis(fig[1,3], xlabel="x", ylabel="y")

scatter!(ax1, t, x)
scatter!(ax2, t, y)
scatter!(ax3, x, y)

fig

pst = posterior(gpxt, y)

ym, yvar = mean_and_var(pst, RowVecs(hcat(x, t)))

ord2 = sortperm(t)
lines!(ax2, t[ord2], ym[ord2])
band!(ax2, t[ord2], ym[ord2] .- yvar[ord2], ym[ord2] .+ yvar[ord2]; color = ("#E69F00", 0.2))

ord3 = sortperm(x)
lines!(ax3, x[ord3], ym[ord3])
band!(ax3, x[ord3], ym[ord3] .- yvar[ord3], ym[ord3] .+ yvar[ord3]; color = ("#E69F00", 0.2))

fig

x_test = range(-5, 5, length=100)
t_test = range(0, 10, length=100)
x_t_grid = ColVecs(reduce(hcat, collect.(Iterators.product(x_test, t_test))))

ym, yvar = mean_and_var(pst, x_t_grid)
contourf(x_test, t_test, reshape(ym, 100, 100)) # Or another 2d contour you'd like