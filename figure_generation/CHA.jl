using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))



# %% sensitive dependence
using DynamicalSystems, PyPlot, Random
Random.seed!(5)
close("all")
fig = figure(figsize=(2figx/3, figx/2))
ax = gca()
u0 = [10,10.0,10]
lo = Systems.lorenz([10,10.0,10])

for i in 1:3
    u = u0 .+ i*1e-3
    tr = trajectory(lo, 15, u)
    plot(0:0.01:15, tr[:, 1],  c = COLORS[[2,3,1][i]])
end

xlabel("\$t\$",labelpad = -10)
yticks(-15:10:15)
xticks(0:5:15)
ylabel("\$x\$",labelpad = -15)
tight_layout()
subplots_adjust(bottom = 0.15, left = 0.12)
fsave(joinpath(figdir, "sensitive"), fig)

# %% stretching and folding for logistic
lo = Systems.logistic()
X = trajectory(lo, 1000, Ttr = 100)
# close("all")
fig = figure()
# ax = subplot2grid((1,3), (0,0))
ax = subplot(121)
ax.scatter(X[1:end-1], X[2:end], color = COLORS[1], s = 5, label = "\$x_{i+2}\$")
ax.set_xlabel("\$x_{i}\$")
ax.set_yticks([0, 0.5, 1])
ax.scatter(X[1:end-2], X[3:end], color = COLORS[2], s = 5, label = "\$x_{i+2}\$")
ax.legend(markerscale = 5)

using3D()
# ax2 = subplot2grid((1,3), (0,1), projection = "3d", colspan = 2)
ax2 = subplot(122, projection = "3d")
ax2.scatter(X[1:end-2], X[2:end-1], X[3:end], color = COLORS[1], s = 5)

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_zticklabels([])
ax2.set_xlabel("\$x_{i}\$", labelpad = -10)
ax2.set_ylabel("\$x_{i+1}\$", labelpad = -10)
ax2.set_zlabel("\$x_{i+2}\$", labelpad = -10)

# fold arrow
s = (0.15, 0.5, 1.2)
e = (0.85, 0.5, 1.2)
ax2.quiver3D(e..., (s .- e)..., color = COLORS[3])
ax2.text3D(0.5, 0.1, 1.2, "Fold", color = COLORS[3])

# stretch arrow
s = (0.75, 0.35, 0.9)
e = (0.95, 0.15, 0.1)
ax2.quiver3D(e..., (s .- e)..., color = COLORS[5])

s = (0.55, 0.45, 0.1)
e = (0.75, 0.35, 0.9)
ax2.quiver3D(e..., (s .- e)..., color = COLORS[5])
ax2.text3D(1.1, 0.45, 0.1, "Stretch", color = COLORS[5])


tight_layout()
fig.subplots_adjust(top = 0.95, wspace = 0.1)
add_identifiers!(fig)
savefig(joinpath(@__DIR__, "stretching.png"))


# %% Folding in henon map
a = 1.4; b = 0.3
he = Systems.henon(;a=a, b=b)
close("all")
fig = figure()

x0, y0 = columns(trajectory(he, 5000, Ttr = 100))
x1 = x0
y1 = @. 1 - a*x0^2 + y0

x2 = b .* x1
y2 = y1

y3 = x2
x3 = y2

ax1 = subplot(141)
ax1.scatter(x0, y0, color = COLORS[1], s = 1, label = "initial")
ax1.text(0.1, 0.05, "initial", transform=ax1.transAxes)

αprev = 0.75
ax2 = subplot(142)
ax2.scatter(x0, y0, color = COLORS[1], s = 0.25, alpha = αprev)
ax2.scatter(x1, y1, color = COLORS[2], s = 1, label = "stretch/bend")
ax2.text(0.1, 0.05, "bend", transform=ax2.transAxes)
ax2.arrow(0,-0.19, 0, 0.8, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[6], length_includes_head=true)
ax2.arrow(-1.5, 0.0,  0, -0.5, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[6], length_includes_head=true)
ax2.arrow(1.5, 0.0,  0, -0.5, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[6], length_includes_head=true)

ax3 = subplot(143)
ax3.scatter(x1, y1, color = COLORS[2], s = 0.25, alpha = αprev)
ax3.scatter(x2, y2, color = COLORS[3], s = 1)
ax3.text(0.1, 0.05, "squeeze", transform=ax3.transAxes)
ax3.arrow(-1.5, 0.0,  1, 0, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[6], length_includes_head=true)
ax3.arrow(1.5, 0.0,  -1, 0, head_width=0.1, head_length=0.1, linewidth=4, color=COLORS[6], length_includes_head=true)


ax4 = subplot(144)
ax4.scatter(x2, y2, color = COLORS[3], s = 0.25, alpha = αprev)
ax4.scatter(x3, y3, color = COLORS[4], s = 1, label = "reflect")
ax4.text(0.1, 0.05, "reflect", transform=ax4.transAxes)

ax4.plot([-2,2], [-2,2], ls = "dashed", lw = 2, color = COLORS[6])

# ticks and stuff
ax1.set_ylabel("\$y\$")
for ax in (ax1, ax2, ax3, ax4)
    ax.set_xlim(-2, 2)
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks(-1.5:0.5:1.5)
    ax.set_xticks(-1.5:1:1.5)
    ax.grid("on")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_xlabel("\$x\$", labelpad = -10)
    # ax.legend(markerscale = 4, scatterpoints=1)
end
tight_layout()
subplots_adjust(left = 0.05, bottom = 0.1, wspace = 0.1)
add_identifiers!(fig; xloc = 0.01)
# savefig(joinpath(figdir, "stretchinghenon"))

# %% definition of lyapunov
using LinearAlgebra
he = Systems.henon()
X₁ = trajectory(he, 50)
u₂ = get_state(he) .+ 1e-6
X₂ = trajectory(he, 50, u₂)
δ  = norm.(X₂.data .- X₁.data)
λ = lyapunov(he, 10000)

close("all")
figure(figsize=(figx/2,figx/4))
ax = gca()
ax.plot(0:50, log.(δ), c = COLORS[1], label ="\$\\ln(\\delta(n)))\$")
ax.set_yticks(-12:4:0)
ax.set_xlabel("\$n\$", labelpad=-10)
# Lyapunov
ax.plot([0, 25] .+ 5, λ .* [0, 25] .- 13, color = COLORS[2])
ax.text(20, -9, "\$\\lambda\$=$(round(λ;digits=2))", color = COLORS[2])
ax.legend()
xticks(0:15:50)
tight_layout()
subplots_adjust(bottom = 0.2, left = 0.12)
savefig(joinpath(figdir, "lyapunov"))

# %% gali map
using DynamicalSystems, PyPlot
figure()
hh = Systems.henonheiles()
ics = Systems.henonheiles_ics(0.13, 10)
cmap = matplotlib.cm.get_cmap("viridis")
for ic in ics
    psos = poincaresos(hh, (1, 0.0), 10000; u0 = ic)
    GALI, t = gali(hh, 1000, 4; u0 = ic)
    v = clamp(t[end]/1000, 0, 1)
    scatter(psos[:, 2], psos[:, 4], color = cmap(v), s = 2)
end

cb = colorbar()
xticks()
xlabel("\$y\$", labelpad = -10)
yticks(-0.5:0.3:0.5)
ylabel("\$p_y\$", labelpad = -10)
cb.set_ticks([0, 1])
cb.set_label("regularity")
tight_layout()
subplots_adjust(bottom = 0.15, top = 0.95, left = 0.1, wspace = 0.01)
savefig(joinpath(figdir, "gali.png"))
