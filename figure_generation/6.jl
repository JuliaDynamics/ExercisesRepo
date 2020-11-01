using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# %% showing embedding
lo = Systems.lorenz([0, 10, 0.0])

tr = trajectory(lo, 1000; Ttr=100)
x, y, z = columns(tr)
w = @. x
τ = estimate_delay(w, "mi_min")
ψ = @. z - y

# figure for Cao's method
fig = figure(figsize=(figx/2.5, figx/2.5));
for (s, l) in zip((w, ψ), ("\$w=x\$", "\$w=z-y\$"))
    τ = estimate_delay(s, "mi_min")
    Ds = estimate_dimension(s, τ, 1:6, "afnn")
    plot(2:7, Ds, marker = "o", label = l)
end
xticks(2:7)
xlabel("\$d\$")
ylabel("\$E_{d}\$")
legend(title="measurement")
fig.tight_layout(pad = 0.1)
# fsave(joinpath(figdir, "caodemonstration"), fig)

# boxes = 10 .^ (-3:0.2:1)
# D_A = generalized_dim(0.0, tr, boxes)
# D_R = generalized_dim(0.0, R, boxes)
# @show (D_A, D_R)

# %%
γ = 2
R = reconstruct(w, γ, τ)

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
γ = 3
τ1 = 2
τ2 = τ
τ3 = 2τ2
θ=-30
e=10

fig = figure()
for (i, τ) in enumerate((0, τ1, τ2, τ3))
    R = i == 1 ? tr : reconstruct(w, γ, τ)
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
