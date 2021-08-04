using AbstractGPs
using KernelFunctions
using GaPLAC
using CairoMakie
using AbstractGPs
using KernelFunctions
using AbstractGPsMakie
using Distances

using CSV
using DataFrames

struct CategoricalKernel <: KernelFunctions.SimpleKernel end
KernelFunctions.kappa(::CategoricalKernel, d2::Real) = d2 > 0 ? 0.0 : 1.0
KernelFunctions.metric(::CategoricalKernel) = Euclidean()

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

##

k_t = with_lengthscale(SqExponentialKernel(), 1)
k_sub = CategoricalKernel()
k_nutrient = LinearKernel()

k1 = (k_t ⊗ k_sub) ∘ SelectTransform([1,2]) + k_nutrient ∘ SelectTransform([3]) # Collect all the kernels to make them act dimension wise
k2 = (k_t ⊗ k_sub) ∘ SelectTransform([1,2]) # kernel without diet variable

# k1 = ((k_t ∘ SelectTransform([1])) ⊗ (k_sub ∘ SelectTransform([2]))) + k_nutrient ∘ SelectTransform([3]) # Collect all the kernels to make them act dimension wise
# k2 = (k_t ⊗ k_sub) ∘ SelectTransform([1,2]) # kernel without diet variable

# Here we create a the prior based on the kernel and the data
priorgp1 = AbstractGPs.FiniteGP(GP(k1), hcat(data.DateOffset, data.PersonID, data.nutrient), 1, obsdim = 1)
priorgp2 = AbstractGPs.FiniteGP(GP(k2), hcat(data.DateOffset, data.PersonID), 1, obsdim = 1)


# We finally compute the posterior given the y observations

posteriorgp1 = posterior(priorgp1, data.bug)
posteriorgp2 = posterior(priorgp2, data.bug)

p1 = logpdf(priorgp1, data.bug)
p2 = logpdf(priorgp2, data.bug)

p1 - p2


##

person_refs = rand(unique(data.PersonID), 10)
dfs = [filter(:PersonID=> ==(p), data) for p in person_refs]
ids = findall(df-> nrow(df) >=3, dfs)
person_refs = person_refs[ids]
dfs = dfs[ids]

dates = range((extrema(data.DateOffset) .+ (-5, 5))..., length = 40)

fig = Figure(resolution=(800,1500));

for (i,p) in enumerate(person_refs)
    testgrid1 = GaPLAC._make_test_rows(dates, [p], dfs[i].nutrient)
    testgrid2 = GaPLAC._make_test_rows(dates, [p])

    ax1 = Axis(fig[i, 1], title = "PersonID $p, with diet", xlabel="Date (offset)", ylabel="bug")
    ax2 = Axis(fig[i, 2], title = "PersonID $p, no diet", xlabel="Date (offset)", ylabel="bug")
    
    ym1, yvar1 = mean_and_var(posteriorgp1, testgrid1)
    ym2, yvar2 = mean_and_var(posteriorgp2, testgrid2)

    n = nrow(dfs[i])
    l = length(dates)

    for d in 1:n
        idx = (1 + (d-1) * l):d*l
        scatter!(ax1, dates, ym1[idx])
        band!(ax1, dates, ym1[idx] .- yvar1[idx], ym1[idx] .+ yvar1[idx], color = (:gray, 0.2))
    end


    scatter!(ax2, dates, ym2, color=:dodgerblue)
    scatter!(ax2, dfs[i].DateOffset, dfs[i].bug, color=:red)
    band!(ax2, dates, ym2 .- yvar2, ym2 .+ yvar2, color=(:dodgerblue, 0.3))
end

fig

##

data = CSV.read("test/testin/input_pair_109.tsv", DataFrame)
filter!(:Date=> !ismissing, data)
data.Date = disallowmissing(data.Date)
gdf = groupby(data, :PersonID)

# Give each subject a random offset from 0
data = transform(gdf, :Date=> (d-> d .- minimum(d) .+ rand()) => :DateOffset)

k_t = with_lengthscale(SqExponentialKernel(), 0.5)
k_sub = CategoricalKernel()
k_nutrient = LinearKernel()

k1 = (k_t ⊗ k_sub) ∘ SelectTransform([1,2]) + k_nutrient ∘ SelectTransform([3]) # Collect all the kernels to make them act dimension wise
k2 = (k_t ⊗ k_sub) ∘ SelectTransform([1,2]) # kernel without diet variable

# Here we create a the prior based on the kernel and the data
priorgp1 = AbstractGPs.FiniteGP(GP(k1), hcat(data.DateOffset, data.PersonID, data.nutrient), 1, obsdim = 1)
priorgp2 = AbstractGPs.FiniteGP(GP(k2), hcat(data.DateOffset, data.PersonID), 1, obsdim = 1)


# We finally compute the posterior given the y observations

posteriorgp1 = posterior(priorgp1, data.bug)
posteriorgp2 = posterior(priorgp2, data.bug)

p1 = logpdf(priorgp1, data.bug)
p2 = logpdf(priorgp2, data.bug)

p1 - p2

##

using StatsBase


mc = mean_and_cov(posteriorgp, RowVecs(hcat(data.PersonID, data.DateOffset, data.nutrient)))

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


##

n = 60

t₀ = collect(range(1, 20, length=n))
f(x) = x + randn() / 10
x₀ = f.(repeat(range(-2, 2, length = n ÷ 3), 3))

y₁ = sin.(t₀) .+ randn(n) ./ 4
y₂ = sin.(t₀ .+ π) .+ randn(n) ./ 4
y₃ = sin.(t₀) .+ 1 .+ randn(n) ./4

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[1,3])
ax4 = Axis(fig[2,1:3])

plot!(ax1, t₀, y₁, color=:blue)
plot!(ax2, t₀, y₂, color=:orange)
plot!(ax3, t₀, y₃, color=:purple)

plot!(ax1, t₀, x₀, color=:grey)
plot!(ax2, t₀, x₀, color=:grey)
plot!(ax3, t₀, x₀, color=:grey)

plot!(ax4, t₀, y₁ .+ x₀, color=:blue)
plot!(ax4, t₀, y₂ .+ x₀, color=:orange)
plot!(ax4, t₀, y₃ .+ x₀, color=:purple)

fig

##

k1 = SqExponentialKernel()
k2 = CategoricalKernel()
k3 = LinearKernel()

gp = GP((k1  ⊗ k2) ∘ SelectTransform([1,2]) + (k3 ∘ SelectTransform([3])))
sub = repeat([1,2,3], inner=n)
t = repeat(t₀, outer=3)
x = repeat(x₀, outer=3)


gptsx = AbstractGPs.FiniteGP(gp, RowVecs(hcat(t, sub, x)), 0.1)

ptsx = posterior(gptsx, [y₁; y₂; y₃])

refsub = 1.0
refx = 0.0
tx_test1 = range(0, 21, length=100)
tx_grid1 = RowVecs(reduce(vcat, hcat.(tx_test1, refsub, refx)))

ym1, yvar1 = mean_and_var(ptsx, tx_grid1)
lines!(ax1, tx_test1, ym1, color=:blue)
band!(ax1, tx_test1, ym1 .+ yvar1, ym1 .- yvar1, color=(:blue, 0.3))
lines!(ax4, tx_test1, ym1, color=:blue)

fig

refsub = 2.0
refx = 0.0
tx_test2 = range(0, 21, length=100)
tx_grid2 = RowVecs(reduce(vcat, hcat.(tx_test2, refsub, refx)))

ym2, yvar2 = mean_and_var(ptsx, tx_grid2)
lines!(ax2, tx_test2, ym2, color=:orange)
band!(ax2, tx_test2, ym2 .+ yvar2, ym2 .- yvar2, color=(:orange, 0.3))
lines!(ax4, tx_test2, ym2, color=:orange)

refsub = 3.0
refx = 0.0
tx_test3 = range(0, 21, length=100)
tx_grid3 = RowVecs(reduce(vcat, hcat.(tx_test3, refsub, refx)))

ym3, yvar3 = mean_and_var(ptsx, tx_grid3)
lines!(ax3, tx_test3, ym3, color=:purple)
band!(ax3, tx_test3, ym3 .+ yvar3, ym3 .- yvar3, color=(:purple, 0.3))
lines!(ax4, tx_test3, ym3, color=:purple)

fig

##

xax = Axis(fig[3,1:3])

for sub in [1., 2., 3.]
    for t in range(0,21, length = 4)
        x = range(-3, 3, length=100)
        st_grid = RowVecs(reduce(vcat, hcat.(t, sub, x)))

        ym, yvar = mean_and_var(ptsx, st_grid)

        lines!(xax, x, ym, color=(:black, t / 40))
        band!(xax, x, ym .+ yvar, ym .- yvar, color=(:black, t / 100))
    end
end

fig

##

using AbstractGPs
using KernelFunctions
using StatsBase


k1 = SqExponentialKernel()
k2 = LinearKernel()

n = 30
t1 = collect(range(1, 20, length=n))
y = randn(n)

f(x) = x + randn() / 10
x = f.(repeat(range(-2, 2, length = n ÷ 3), 3))

gp = GP((k1 ∘ SelectTransform([1])) + (k2 ∘ SelectTransform([2])))
fgp = AbstractGPs.FiniteGP(gp, RowVecs(hcat(t1, x)), 0.1)

gppost = posterior(fgp, y)

testt = range(20,40, length=100)
refx = 0
testgrid = ColVecs(reduce(hcat, vcat.(testt, refx)))
ym1, yvar1 = mean_and_var(gppost, testgrid)

##

using Distributions

x1 = 1:20
x2 = randn(20)

k1 = SqExponentialKernel()
k2 = LinearKernel()

K = (k1 ⊗ k2) ∘ SelectTransform([1,2])
pr1 = AbstractGPs.FiniteGP(GP(K), hcat(x1,x2), obsdim=1)
rand(pr1)

K = (k1 ∘ SelectTransform(1)) ⊗ (k2 ∘ SelectTransform(2))
pr2 = AbstractGPs.FiniteGP(GP(K), hcat(x1,x2), obsdim=1)
rand(pr2)

