using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

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

# fsave(joinpath(figdir, "permentropy"), fig)


# %% Koch snowflake
linepoints = [[0.0; 0.0], [1.0; 0.0]]
flakepoints = [[0.0; 0.0], [0.5; sqrt(3)/2], [1; 0.0], [0.0; 0.0]]
function koch(points, maxk, α = sqrt(3)/2)
  Q = [0 -1; 1 0]
  for k = 1:maxk
    n = length(points)
    new_points = Vector{Float64}[]
    for i = 1:n-1
      p1, p2 = points[i], points[i+1]
      v = (p2 - p1) / 3
      q1 = p1 + v
      q2 = p1 + 1.5v + α * Q * v
      q3 = q1 + v
      append!(new_points, [p1, q1, q2, q3])
    end
    push!(new_points, points[end])
    points = new_points
  end
  return points
end

# Plot construction
fig, axes = subplots(2,3, figsize = (2figx/3,figx/2 - 2))
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
    botrow[i].add_artist(plt.Circle((0.5,0.28900), 0.5, lw=2.0, color="C1", alpha = 0.75, fill=false))
    botrow[i].axis("off")
end
for i in 1:2
    origin = botrow[i]
    zoomin = botrow[i+1]
    zbox = allpoints[i+1]
    axis_zoomin!(zoomin, origin, zbox, zbox, "C$(i+1)", 0.4)
end

fig.tight_layout()
fig.subplots_adjust(wspace = 0.01, hspace = 0.01, bottom = 0.05, left = 0.05)
fsave(joinpath(figdir, "koch"), fig)

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
# fig.savefig(joinpath(figdir, "henon_zoom.png"))


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


# %% add
fig = figure(;figsize = (figx/2, figy))
ax = gca()
he = Systems.henon()
N = 10000
tr = trajectory(he, N, Ttr = 100)
es = [0.2, 0.05, 0.01]
ax.grid(false)
mini = minima(tr)

for (i, e) ∈ enumerate(es)
    for x in mini[1]:e:1.4
        ax.axvline(x; color = "k", lw = 1.0, alpha = 1/(3^i))
    end
    for y in mini[2]:e:0.4
        ax.axhline(y; color = "k", lw = 1.0, alpha = 1/(3^i))
    end
    p, bins = binhist(e, tr)
    for b in bins
        r = matplotlib.patches.Rectangle(b, e, e; alpha = 0.8, color = "C$i")
        ax.add_artist(r)
    end
end
ax.plot(tr[:, 1], tr[:, 2], ls = "None", marker = ".", color = COLORS[1], ms = 1, zorder = 99)
ax.set_yticks([])
ax.set_ylabel("\$y\$")
ax.set_xticks([])
ax.set_xlabel("\$x\$")
fig.subplots_adjust(left=0.05, bottom = 0.1, right = 0.98, top = 0.98)
wsave(plotsdir("henon_gridding"), fig)
