using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# %% Koch snowflake
linepoints = SVector{2}.([[0.0; 0.0], [1.0; 0.0]])
flakepoints = SVector{2}.([[0.0; 0.0], [0.5; sqrt(3)/2], [1; 0.0], [0.0; 0.0]])
function koch(points, maxk, α = sqrt(3)/2)
  Q = @SMatrix [0 -1; 1 0]
  for k = 1:maxk
    n = length(points)
    new_points = eltype(points)[]
    for i = 1:n-1
      p1, p2 = points[i], points[i+1]
      v = (p2 - p1) / 3
      q1 = p1 + v
      q2 = p1 + 1.5v + α * Q * v
      q3 = q1 + v
      push!(new_points, p1, q1, q2, q3)
    end
    push!(new_points, points[end])
    points = new_points
  end
  return points
end

# Plot construction
fig, axes = subplots(2,3, figsize = (figx, 2figy))
toprow = axes[1:2:5]
botrow = axes[2:2:6]

for (o, ax) in enumerate(toprow)
    ax.clear()
    ax.set_aspect("equal")
    k1 = koch(flakepoints, o-1)
    k2 = koch(flakepoints, o)

    ax.plot([a[1] for a in k1], 0.5 .+ [a[2] for a in k1], COLORS[o])
    ax.plot([a[1] for a in k2], 0.5 .+ [a[2] for a in k2], COLORS[o+1], label="Step $(o)", alpha = 0.75)
    # ax.legend(loc = "upper left", framealpha = 1.0, ncol = 2, handlelength = 1)
    ax.text(-0.1, 0.88, "step $(o)", transform = ax.transAxes, size=18)
    # ax.set_ylim(ys...)
    # ax.set_xlim(xs...)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis("off")
end

# Plot versus circle
points1 = [[0.5; 0.865], [2/3, 0.577361]]
points2 = [[0.537233; 0.801228], [0.555887, 0.768912]]
allpoints = [flakepoints, points1, points2]

for i in 1:3
    botrow[i].clear()
    largekoch = koch(allpoints[i], 6)
    kx, ky = [a[1] for a in largekoch], [a[2] for a in largekoch]
    botrow[i].plot(kx, ky, lw = 0.5)
    botrow[i].set_aspect("equal")
    botrow[i].add_artist(plt.Circle(
        (0.5,0.28900), 0.5, lw=2.0, color="C1", alpha = 0.75, fill=false
    ))
    botrow[i].axis("off")
end
for i in 1:2
    origin = botrow[i]
    zoomin = botrow[i+1]
    zbox = allpoints[i+1]
    axis_zoomin!(zoomin, origin, zbox, zbox, "C$(i+1)"; α = 0.5)
end

fig.tight_layout(pad = 0.1)
fig.subplots_adjust(wspace = 0.01, hspace = 0.01, bottom = 0.01, left = 0.01)
# save(plotsdir("koch"), fig)

# %% Henon fractal zoom
fig, axes = subplots(1,3, figsize = (figx,figx/4))
botrow = axes

he = Systems.henon()
N = 10000
tr = trajectory(he, N, Ttr = 100)
integ = integrator(he)
data1 = Dataset([integ.u])
data2 = Dataset([integ.u])
n, m = 1, 1

zbox1 = ((0.93, 0.05), (1.17, 0.13))
zbox2 = ((1.1275, 0.078), (1.1425, 0.085))
while n < N
    step!(integ)
    if zbox1[1][1] < integ.u[1] < zbox1[2][1] && zbox1[1][2] < integ.u[2] < zbox1[2][2] && m < N
        push!(data1, integ.u)
        global m+=1
    end
    if zbox2[1][1] < integ.u[1] < zbox2[2][1] && zbox2[1][2] < integ.u[2] < zbox2[2][2]
        push!(data2, integ.u)
        global n+=1
    end
end

for ax in botrow[1:2]
    ax.clear()
    ax.plot(tr[:, 1], tr[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1)
end
botrow[1].set_xticks([])
botrow[1].set_yticks([])
botrow[2].clear()
botrow[2].plot(data1[:, 1], data1[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1)
botrow[2].axis("off")
botrow[3].clear()
botrow[3].plot(data2[:, 1], data2[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1)
botrow[3].axis("off")

for (i, zbox) in enumerate((zbox1, zbox2))
    axis_zoomin!(botrow[i+1], botrow[i], zbox, zbox, COLORS[i+1])
    botrow[i+1].set_xlim(zbox[1][1], zbox[2][1])
    botrow[i+1].set_ylim(zbox[1][2], zbox[2][2])
end

fig.tight_layout()
fig.subplots_adjust(top = 0.98, bottom = 0.05, right = 0.95, left = 0.05, wspace=0.1, hspace = 0.1)
# wsave(plotsdir("henon_zoom"), fig)


# %% dimension: Entropy and correlation
he = Systems.henon()
N, α = 10000, 2
X = trajectory(he, N, Ttr = 10)
ε = 10 .^ (-5.0:0.1:1)
H = genentropy(α, ε, X)
C = correlationsum(X, ε)

fig, axs = subplots(2, 1; figsize = (figx/2, figx/2), sharex = true)

axs[1].plot(log.(ε), -H, COLORS[1])

x = -5:0
D = 1.22


axs[2].plot(log.(ε), log.(C), COLORS[1])


for (ax, x) in zip(axs, (-5:0, -10:0))
    Dl = D .* x
    Dl = Dl .- 2
    ax.plot(x, Dl, c = COLORS[3])
    ax.axvspan(x[1], x[end], color = COLORS[2], alpha = 0.25)
    ax.grid("on")
    ax.text(x[2] + 1.5, Dl[2], "\$D_2\$ = $D", color = COLORS[3], size = 30)

end

axs[2].set_xticks(-10:2:2)
axs[1].tick_params(labelbottom=false)
axs[2].set_xlabel("\$\\log(\\varepsilon)\$")
axs[1].set_ylabel("\$-H_$α\$", labelpad = -20)
axs[2].set_ylabel("\$\\log(C)\$", labelpad = -20)
axs[1].set_yticks(-10:2:0)
axs[1].set_yticklabels(["-10", "", "", "", "", "0"])
axs[2].set_yticks(0:-2:-12)
axs[2].set_yticklabels(["0", "", "", "",  "","", "-12"])
fig.tight_layout()
fig.subplots_adjust(left=0.14, bottom = 0.15, hspace = 0.1)
# fsave(joinpath(figdir, "dimension"), fig)


# %% Henon gridding for dimension explanation
fig = figure()
ax = subplot(1,2,1)
ax2 = subplot(1,2,2)
he = Systems.henon()
N = 10000
tr = trajectory(he, N, Ttr = 100)
es = [0.2, 0.05, 0.01]
ns = zero(es)
ax.grid(false)
mini = minima(tr)

for (i, e) ∈ enumerate(es)
    if i ≤ 2
        for x in mini[1]:e:1.4
            ax.axvline(x; color = "k", lw = 0.5, alpha = 1/(3^i))
        end
        for y in mini[2]:e:0.4
            ax.axhline(y; color = "k", lw = 0.5, alpha = 1/(3^i))
        end
    end
    p, bins = binhist(tr, e)
    for b in bins
        r = matplotlib.patches.Rectangle(b, e, e; color = "C$i", ec = "k", lw = 1/i)
        ax.add_artist(r)
    end
    ns[i] = length(probabilities(tr, e))
end

ax.plot(tr[:, 1], tr[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1, zorder = 99)
ax.set_yticks([])
ax.set_ylabel("\$y\$")
ax.set_xticks([])
ax.set_xlabel("\$x\$")

ax2.scatter(log.(1 ./ es), log.(ns); c = ["C$i" for i in 1:length(es)], s = 200, zorder = 99)
s = linreg(log.(1 ./ es), log.(ns))[2]
ax2.plot(log.(1 ./ es), log.(1 ./ es) .* s .+ 2, color = "C0")
ax2.text(3.5, 6, "\$\\Delta = $(rdspl(s))\$"; color = "C0")
ax2.set_xlabel("\$\\log ( 1/\\varepsilon)\$"; labelpad = -15)
ax2.set_ylabel("\$\\log ( M)\$")
ax2.set_xticks([2, 3, 4])
ax2.set_xticklabels([2, "", 4])
fig.tight_layout(pad = 0.25)
fig.subplots_adjust(wspace = 0.2)
wsave(plotsdir("henon_gridding"), fig)

# %% Noise radius illustration using e.g. Poincare section
using DynamicalSystems, PyPlot
using Printf

ro = Systems.roessler([0.1,0.2,0.1])
fig, axs = subplots(2, 1; figsize = (0.4*figx, 1.5figy))
plane = (2, 0.0)

err = (1e-4, 1e-12)
εs = 10 .^ (-5:0.5:1)

for (i, e) ∈ enumerate(err)
    p = poincaresos(
        ro, plane, 10000; Ttr = 100.0,
        rootkw = (xrtol = e, atol = e)
    )

    axs[1].scatter(p[:, 1], p[:, 3]; s = 20, alpha = 0.75/i^2,
    label = "tol = \$10^{$(round(Int, log10(e)))}\$")
    Cs = correlationsum(p, εs)
    y = log.(Cs ./ maximum(Cs))
    axs[2].plot(log.(εs), y)

    # add noise line
    if i == 1
        z = log.(εs)[end-3]
        axs[2].plot([z, z], [-14, y[end-3]]; color = "C2", lw = 2, ls = ":")
        axs[2].text(z+0.2, -9, "\$\\sigma\$"; color = "C2")
    end
end

axs[1].set_xlabel("\$x\$";labelpad = -15)
axs[1].set_ylabel("\$z\$")
axs[1].legend(markerscale = 5, fontsize=26)
axs[1].set_yticks(0:6:18)
axs[2].plot([-9, 0], 0.9 .* [-9, 0]; ls = "--", color = "C1")
axs[2].text(-7, -4, "0.9"; color = "C1")
axs[2].plot([-6, -3], 2.0 .* [-6, -3] .- 2; ls = "--", color = "C0")
axs[2].text(-3, -11, "2"; color = "C0")
axs[2].set_ylabel("\$\\log(C)\$"; labelpad = -10)
axs[2].set_xlabel("\$\\log(\\varepsilon)\$"; labelpad = -12)
axs[2].set_xticks(-12:5:2)
axs[2].set_yticks(-14:5:2)

fig.tight_layout(pad = 0.5)
add_identifiers!(fig)
wsave(plotsdir("fractaldim_noise"), fig)


# %% Magnetic pendulum
using LinearAlgebra

α=0.2; ω=1.0; d=0.3
gx = gy = range(-5, 5; length = 1500)
config = @strdict(α, ω, d, gx, gy)

function magnetic_basins(config)
    @unpack α, ω, d, gx, gy = config
    ma = Systems.magnetic_pendulum(; α, ω, d)
    @time basins, attractors = basins_general(gx, gy, ma; idxs = 1:2)
    return @strdict(gx, gy, basins, attractors)
end

# produce high quality
produce_or_load(datadir(), config, magnetic_basins; prefix = "magnetic_basins")

# produce zoomed version
gx = range(1.80, 1.95; length = 1000)
gy = range(0, 0.12; length = 1000)
config = @strdict(α, ω, d, gx, gy)
produce_or_load(datadir(), config, magnetic_basins; prefix = "magnetic_basins_zoomed")

# Plot this
gx, gy, basins, attractors = @unpack wload(datadir(savename("magnetic_basins", config)))

LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])

fig = figure(figsize=(figx/2, figy))
pcolormesh(gx, gy, basins'; cmap, shading = "gouraud")
gca().set_aspect("equal")
xticks([-5, 5])
yticks([-5, 5])
xlabel("\$x\$", labelpad=-30)
ylabel("\$y\$", labelpad=-30)
for m in attractors
    m1, m2 = columns(m)
    scatter(m1, m2; color = "white", edgecolor = "black", zorder = 99, s = 100)
end
tight_layout(pad=0.3)
subplots_adjust(bottom = 0.1, top = 0.95)
wsave(plotsdir("magneticpendulum"), fig)

# Plot zoomed version
gx, gy, basins, attractors = @unpack wload(datadir(savename("magnetic_basins_zoomed", config)))
LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])

fig = figure(figsize=(figx/2, figx/2))
pcolormesh(gx, gy, basins'; cmap, shading = "gouraud")
gca().set_aspect("equal")

xticks([1.8, 1.95], size = 20)
yticks([0, 0.12], size = 20)

fig.savefig(joinpath(figdir, "magneticpendulum_zoom.png"))

# xticks([-5, 5])
# yticks([-5, 5])
