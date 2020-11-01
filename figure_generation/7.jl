using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# %% mutual information
using DynamicalSystems, PyPlot
using InformationMeasures, Random
MI = get_mutual_information #shortcut
lo, N = Systems.logistic(), 100
x = trajectory(lo, N-1)
y = trajectory(lo, N-1; Ttr = 1)
x .+= randn(N)/25; y .+= randn(N)/25
m = MI(x, y)
null = [MI(shuffle!(x), shuffle!(y))
        for _ in 1:10000]

# plot stuff
fig = figure(figsize = (0.9figx/2, 0.8figx/2))
using Statistics
μ, σ = mean(null), std(null)
plt.hist(null, 50, label = "null pdf")
axvline(μ, color = "C1", label = "\$\\mu\$", ls = "dashed")
axvline(μ-3σ, color = "C2", label = "\$\\mu \\pm 3\\sigma\$", ls = "dashed")
axvline(μ+3σ, color = "C2", ls = "dashed")
axvline(m, color = "C3", label = "\$m\$")

legend(framealpha = 1.0, loc = "upper right")
yticks([])
xlabel("mutual information (bits)")
tight_layout()
# fsave(joinpath(figdir, "mutualinfo"), fig)

# %% transfer entropy
using DynamicalSystems, CausalityTools, PyPlot, Random

function ulam(dx, x, p, t)
    f(x) = 2 - x^2; ε = p[1]; N = length(x)
    for i in 1:N
        dx[i] = f(ε*x[mod1(i-1, N)] + (1-ε)*x[i])
    end
end
ds = DiscreteDynamicalSystem(ulam, rand(100), [0.04])

methods = [RectangularBinning(r) for r in (0.01, 0.1, 0.4)]
εs = 0.0:0.01:1.0
tes = [zeros(length(εs), 2) for j in 1:length(methods)]

for (i, ε) in enumerate(εs), (j, meth) in enumerate(methods)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, 10000; Ttr = 10000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    tes[j][i, 1] = transferentropy(X1, X2, meth, 1, 1, 1)
    tes[j][i, 2] = transferentropy(X2, X1, meth, 1, 1, 1)
end



fig = figure(figsize = (6figx/10, figx/3))
rs = (0.01, 0.1, 0.4)

for j in 1:length(methods)
    plot(εs, tes[j][:, 1], color = "C$(j-1)", label = "\$r=$(rs[j])\$")
    plot(εs, tes[j][:, 2], color = "C$(j-1)", ls = "dashed")
end
xlabel("coupling strength \$\\epsilon\$")
ylabel("transfer entropy")
legend(loc = "center")
tight_layout()
# fsave(joinpath(figdir, "transfer"), fig)

# %% Surrogate application (From Lancaster)
using DynamicalSystems, TimeseriesSurrogates,
      PyPlot, Statistics, Random

fig = figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,4,3)
ax3 = fig.add_subplot(1,4,4)
axs = [ax1, ax2, ax3]

ro = Systems.roessler(ones(3); a = 0.165, b = 0.2, c = 10.0)

tr = trajectory(ro, 500; dt = 0.1, Ttr = 500)
x = tr[:, 1]
x ./= std(x)
x .+= randn(length(x))*std(x)*0.1
axs[1].plot(x .+ 2,  label = "Rössler", lw = 1.2)

# generate autoregressive process
Random.seed!(77163)
η = randn(5000)
s = ones(5000)
for n in 4:5000
    s[n] = 1.625s[n-1] - 0.284s[n-2] - 0.355s[n-3] + η[n] - 0.96η[n-1]
end
s ./= std(s)

# # equivalent ARFIMA (doesn't work)
# using ARFIMA
# φ = SVector(1.625, -0.284, -0.355)
# θ = SVector(0.96)
# Random.seed!(77163)
# s2 = arma(5000, 1.0, φ, θ)
# s2 ./= std(s2)
# plot(s)
# plot(s2)
# xlim(0, 1000)

extrema(s)
axs[1].plot(s .- 2,label = "ARMA", lw = 1.0)
axs[1].set_xlim(0, 1000)
axs[1].set_ylim(-5, 5)
axs[1].set_yticks(-4:2:4)
# axs[1].set_xticklabels([])
leg = axs[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left",
           ncol=2, mode="expand", borderaxespad=0.)
# leg = axs[1].legend(ncol=1, loc = "center right", fontsize = 24)
for line in leg.get_lines()
    line.set_linewidth(4.0)
end
# axs[1].set_ylabel("timeseries")

# Do the surrogate calculation
εro, εma = std.((x, s))./4
algs = [RandomFourier(), AAFT(), IAAFT()]
names = ["FT", "AAFT", "IAAFT"]
sgx = [surrogenerator(x, m) for m in algs]
sgs = [surrogenerator(s, m) for m in algs]
τx = estimate_delay(x, "ac_zero")
τs = estimate_delay(s, "ac_zero")
Cx = takens_best_estimate(embed(x, 3, τx), εro)
Cs = takens_best_estimate(embed(s, 4, τs), εma)
A = length(algs)

xboxes = []
sboxes = []
for i in 1:A
    sx, ss = sgx[i], sgs[i]
    bx, bs = [], []
    for j in 1:100
        X = embed(sx(), 3, τx)
        S = embed(ss(), 4, τs)
        push!(bx, takens_best_estimate(X, εro))
        push!(bs, takens_best_estimate(S, εma))
    end
    push!(xboxes, bx)
    push!(sboxes, bs)
end

for (j, boxes) in enumerate((xboxes, sboxes))
    c = "C$(j-1)"
    ax = axs[j+1]
    for (i, b) in enumerate(boxes)
        ax.boxplot(b; positions = [i],
        patch_artist=true, boxprops=Dict(:facecolor=>c, :color=>c),
        medianprops=Dict(:color=>"w"), flierprops=Dict(:markeredgecolor=>c))
    end
end

l = [1-0.2, A+0.2]
axs[2].plot(l, fill(Cx, 2), color = "C0", ls = "dashed")
# axs[2].plot(l, fill(2.75, 2), color = "C0", ls = "dashed")
axs[3].plot(l, fill(Cs, 2), color = "C1", ls = "dashed")
# axs[2].set_ylim(2.7, 3)
axs[2].set_yticks(2.6:0.2:3.0)
axs[2].set_title("\$D_C\$, Rössler")
axs[3].set_title("\$D_C\$, ARMA")
axs[3].set_yticks(3.7:0.1:3.8)

for ax in axs[2:3]
    ax.grid(false; axis = "x")
    ax.set_xticks(1:A)
    ax.set_xticklabels(names, rotation = 0, size = 20)
end
fig.subplots_adjust(top = 0.85, left = 0.08, bottom = 0.1, right = 0.98, wspace = 0.3)
# fsave(joinpath(figdir, "surrogates"), fig)



# %% Convergent Cross Mapping
using DynamicalSystems, Neighborhood,
Statistics, LinearAlgebra
function ccm(x, y, d, τ)
    Mx = embed(x, d, τ); theiler = Theiler(0); tree = KDTree(Mx)
    idxs = bulkisearch(tree, Mx.data, NeighborNumber(d+1), theiler)
    ỹ = copy(y)
    for i in 1:length(Mx)
        J = idxs[i]
        xᵢ = Mx[i]; n1 = norm(xᵢ - Mx[J[1]])
        w = [exp(-norm(xᵢ - Mx[j])/n1) for j in J]
        w ./= sum(w)
        ỹ[i] = sum(w[k]*y[j] for (k, j) in enumerate(J))
    end
    return cor(y, ỹ)
end


ds = Systems.lorenz()
tr = trajectory(ds, 5000; dt = 0.1)
x, y = columns(tr)
Ns = 5000:5000:length(tr)
ρs = [ccm(x[1:N], y[1:N], 3, 5) for N in Ns]

ccm(x, y, 3, 5)

plot(Ns, ρs)
