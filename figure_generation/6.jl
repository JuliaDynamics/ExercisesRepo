using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# %% Demonstration of delay embeddings
lo = Systems.lorenz([0, 10, 0.0])
tr = trajectory(lo, 1000; Ttr=10)
x, y, z = columns(tr)
w = x
τ = estimate_delay(w, "mi_min")
R = embed(w, 3, τ)

close("all")
N = 5000
fig = figure()
ax1 = subplot(131, projection = "3d")
ax1.plot3D(x[1:N], y[1:N], z[1:N], lw = 1.0, color = "C1")
ax1.set_title("original set \$A\$")
ax = ax1
for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end
# ax1.text(-20, -20, 0, "\$D_0(A)\$ = $(round(D_A;digits=3))", color = "C1")

ax2 = subplot(132)
wx = 400:1000
ax2.plot(wx, w[wx])
tau = 700:(700+τ)
ax2.plot(tau, fill(10, length(tau)), color = "C3")
ax2.text(700, 7.5, "\$\\tau\$", color = "C3", size = 32)
ax2.set_title("measurement \$w=x\$")
# make fancy axis
ax2.spines["bottom"].set_position("center")
ax2.spines["right"].set_color("none")
ax2.spines["top"].set_color("none")
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel("\$t\$")
ax2.set_ylabel("\$w\$", rotation = 0)
ax2.xaxis.set_label_coords(1.0, 0.5)
ax2.yaxis.set_label_coords(-0.05, 0.9)

ax3 = subplot(133, projection = "3d")
ax3.plot(R[1:N, 1], R[1:N, 2], R[1:N, 3], color = "C2", lw = 1.0)
ax = ax3
for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end
ax3.set_title("reconstruction \$R\$")
# ax3.set_title("reconstruction \$R\$\n\$(\\gamma=$(γ),\\tau=$(τ))\$")
# ax3.text(-5,-25,-10, "\$D_0(R)\$ = $(round(D_R;digits=3))", color = "C2")

fig.tight_layout(pad = 0.1)
# fsave(joinpath(figdir, "delayembedding"), fig)

# %% value of τ demonstration
d = 2
τ1 = 2
τ2 = τ
τ3 = 2τ2
θ=-30
e=10

fig = figure()
for (i, τ) in enumerate((0, τ1, τ2, τ3))
    R = embed(w, d, τ)
    ax = subplot(1, 4, i, projection = "3d")
    ax.plot(R[1:N, 1], R[1:N, 2], R[1:N, 3],
    color = i == 1 ? "C1" : "C2", lw = 1,
    zorder = 1)
    ax.view_init(elev=e, azim=θ)
    ax.set_title(i == 1 ? "original" : "\$R:\\,\\tau=$(τ)\$", pad = -10)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end
    dj = 65
    for j in (456, )
        χ, ψ, ζ = R[j:dj:j+dj, 1], R[j:dj:j+dj, 2], R[j:dj:j+dj, 3]
        ax.scatter3D(χ, ψ, ζ,depthshade = false,c = "C0", zorder = 99, s = 50)
    end
end
# ax = subplot(1, 4, 1, projection = "3d")
# ax.plot3D(x[1:N], y[1:N], z[1:N], lw = 1.0, color = "C1")
# ax.view_init(elev=e, azim=θ)
# for a in (ax.xaxis, ax.yaxis, ax.zaxis); a.set_ticklabels([]); end

fig.tight_layout()
fig.subplots_adjust(wspace = 0.0001, bottom = -0.1, top = 1, left = 0.0, right = 1)
# fsave(joinpath(figdir, "taudemonstration"), fig)


# %% show AC, MI for choosing optimal τ
ττ = 0:90
ac = autocor(w, ττ)
smi = selfmutualinfo(w, ττ)
smi ./= maximum(smi)

fig = figure(figsize = (figx/2, figy))
ax = gca()
ax.plot(ττ, ac; label = "AC")
ax.plot(ττ, smi, label = "SMI")
τ = estimate_delay(w, "mi_min")
ax.scatter(τ, smi[τ]; color = "C2", s = 100, zorder = 99)
ax.set_ylim(0,1.05)
ax.set_xticks(0:30:90)
ax.set_xlabel("\$\\tau\$", labelpad = -20)
ax.legend()
fig.tight_layout(;pad = 0.25)
wsave(plotsdir("mutualinfo_ac_tau"), fig)

# %% Using Cao's method to estimate embedding
lo = Systems.lorenz([0, 10, 0.0])

tr = trajectory(lo, 1000; Ttr=10)
x, y, z = columns(tr)
fig = figure(figsize=(figx/2, figy));
w = x
ψ = z .- y

for (s, l) in zip((w, ψ), ("\$w=x\$", "\$w=z-y\$"))
    τ = estimate_delay(s, "mi_min")
    Ed = delay_afnn(s, τ, 2:7)
    plot(2:7, Ed, marker = "o", label = l)
end

xlabel("\$d\$"; labelpad = -10)
ylabel("\$E_{d}\$")
legend(title="measurement")
fig.tight_layout(pad = 0.4)
fig.subplots_adjust(left = 0.15, bottom = 0.15)
wsave(plotsdir("caodemonstration"), fig)


# %% broomhead king
using DynamicalSystems, PyPlot
ds = Systems.lorenz()
tr = trajectory(ds, 100.0; Ttr=10)
w = tr[:, 1]
w .+= rand(length(w)) #add noise
U, S = broomhead_king(w, 42)
R = reconstruct(w, 2, 17)

close("all")
fig = figure(figsize = (1.2figx/2, figx/2))
subplot(121)
plot(U[:, 1], U[:, 2], lw = 1.0)
title("Broomhead-King")
xticks([]); yticks([])
subplot(122)
plot(R[:, 1], R[:, 2], lw = 1.0, color = "C2")
xticks([]); yticks([])
title("delay embedding")
tight_layout()
subplots_adjust(wspace = 0.05, left = 0.025, right = 0.975)
fsave(joinpath(figdir, "broomhead"), fig)

# %% Comparison of fractal dim of embedding and reconstruction
ds = Systems.lorenz()
A = trajectory(ds, 500; dt = 0.05, Ttr = 1000)
x = A[:, 1]
τ = estimate_delay(x, "mi_min")
d = 3 # embedding dimension = 3, we know it already
R = embed(x, d, τ)
εs = MathConstants.e .^ (-5:0.5:4)
CA = correlationsum(A, εs)
CR = correlationsum(R, εs)

e = log.(εs)
fig = figure(figsize = (figx/2, figy))
PyPlot.plot(e, log.(CA); label = "original", color = "C1")
PyPlot.plot(e, log.(CR); label = "reconstructed", color = "C2")
xlabel("\$\\log(\\varepsilon)\$")
ylabel("\$\\log(C)\$")
PyPlot.legend()
fig.tight_layout(pad = 0.3)


# %% permutation entropy
using Random, Combinatorics, Statistics
Random.seed!(356)
o = 3
N = 8
x = rand(N)

fig = figure(figsize = (figx, figx/3.5))
ax1 = subplot2grid((1, 3), (0, 0))
ax1.plot(1:N, x, marker = "o", color = "C0", mfc = "C1", ms = 15, zorder = 99)
ax1.set_title("timeseries", size = 28)
ax1.spines["bottom"].set_position("center")
ax1.spines["right"].set_color("none")
ax1.spines["top"].set_color("none")
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel("\$t\$")
ax1.set_ylabel("\$x\$", rotation = 0)
ax1.xaxis.set_label_coords(1.0, 0.5)
ax1.yaxis.set_label_coords(-0.05, 0.9)
ax1.set_ylim(-0.1, 1.1)
# ax1.axis("off")

# add pattern indications
for (s, n) in zip(("2", "3", "3"), (2, 4, 6))
    xs = [n, n, n+o-1, n+o-1]
    ym = minimum(x[xs[1]:xs[end]]) - 0.1
    yd = ym - 0.02
    ys = [ym, yd, yd, ym] #.+ 0.05
    ax1.plot(xs, ys, color = "k", lw = 1)
    ax1.text(mean(xs), minimum(ys)-0.02, "#=$s", ha="center", va = "top")
end

# ax1.annotate("4", (2, -0.2), xytext = (4, -0.2), ha="center")

ax2 = subplot2grid((1, 3), (0, 1), colspan = 2)
p = permutations([0, 0.5, 1], o) |> collect
counts = [0, 1, 3, 1, 0, 1]
for (i, a) in enumerate(p)
    ax2.plot((1:o) .+ i*o .+ 1, a, marker = "o",
    color = "C0", mfc = "C2", ms = 10)
    ax2.text(o÷2 + i*o + 1, 1.2, "$i")
    ax2.text(o÷2 + i*o + 1, -0.5, "$(counts[i])")
end
ax2.text(o+1, 1.2, "#", ha = "right")
ax2.text(o+1, 0.5, "pattern", ha = "right")
ax2.text(o+1, -0.5, "count", ha = "right")
ax2.set_ylim(-1, 1.5)
ax2.set_xlim(0, ax2.get_xlim()[2])
ax2.set_title("relative amplitude permutations, \$d=$o\$", size = 28)
ax2.axis("off")

fig.tight_layout()
fig.subplots_adjust(top = 0.9, bottom = 0.02, right = 0.95, left = 0.05, wspace=0.1, hspace = 0.1)

# wsave(plotsdir("permentropy"), fig)
