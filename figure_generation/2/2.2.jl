# Fixed point classification in 2D
using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

alleigs = []
fig, axs = subplots(2,3; figsize = (figx, 2figy))
ax = gca()
xgrid = -5:0.05:5
ygrid = xgrid
ux = zeros(length(xgrid), length(ygrid))
uy = copy(ux)

Ms = [
     [1.0  0.8;
      0.4 1.0],
     [-1.0  -0.8;
      1.8  -1.0],
     [-1.0  0;
      1  1],
     [0 -1;
      1  0],
]

titles = ["repulsive node", "attractive spiral", "hyperbolic/saddle", "center"]

using LinearAlgebra
function stream_eigs!(ax, M, c = "C0")
    for (i, x) in enumerate(xgrid)
        for (j, y) in enumerate(ygrid)
            ux[i, j], uy[i, j] = M * [x, y]
        end
    end

    ax.streamplot(Vector(xgrid), Vector(ygrid), ux', uy';
        linewidth = 1.5, density = 0.5, color = c, arrowsize = 2,
    )
    ev = eigen(M)
    push!(alleigs, ev.values)
    if eltype(ev.values) <: Float64
        e1 = ev.vectors[:, 1]
        e2 = ev.vectors[:, 2]
        for e in (e1, e2)
            ax.plot(e[1] .* 2xgrid, e[2] .* 2ygrid; color = "C0", ls = "dashed")
        end
    end
    ax.set_xlim(xgrid[1], xgrid[end])
    ax.set_ylim(ygrid[1], ygrid[end])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.plot(0, 0; marker = "o", mec = "C0", mew = 1,
        markersize = 12, mfc = "C0", zorder = 99
    )
end

function tovec(e)
    if eltype(e) <: Real
        e1 = e
        e2 = zeros(2)
    else
        e1 = [e[1].re, e[2].re]
        e2 = [e[1].im, e[2].im]
    end
    return e1, e2
end

for (i, M) in enumerate(Ms)
    ax = axs[i]
    stream_eigs!(ax, M, "C$(i)")
    ax.set_title(titles[i]; color = "C$i")
    # Plot eigenvalues
    e1, e2 = tovec(alleigs[i])
    axs[5].scatter(e1, e2; color = "C$i", s = 200, zorder = 99)
end

# Set axis of eigenval plot
axs[5].spines["left"].set_position("center")
axs[5].spines["bottom"].set_position("center")
axs[5].grid(false)
axs[5].spines["right"].set_color("none")
axs[5].spines["top"].set_color("none")
axs[5].set_title("eigenvalues"; color = "k")
# axs[5].set_yticklabels([])
# axs[5].set_xticklabels([])
axs[5].set_xlim(-2, 2)
axs[5].set_ylim(-1.8, 1.8)

# Plot limit cycle
function vanderpoll(u, p, t)
    x, y = u;
    xdot = p*(x - x^3/3 - y)
    ydot = x/p
    return SVector(xdot, ydot)
end

ds = ContinuousDynamicalSystem(vanderpoll, rand(2), 0.5)
tr = trajectory(ds, 1000.0; Ttr = 100.0)

axs[6].plot(columns(tr)...; color = "C0")
for (i, x) in enumerate(xgrid)
    for (j, y) in enumerate(ygrid)
        ux[i, j], uy[i, j] = ds.f(SVector(x, y), ds.p, 0)
    end
end
axs[6].streamplot(Vector(xgrid), Vector(ygrid), ux', uy';
    linewidth = 1.5, density = 0.5, color = "C5", arrowsize = 2
)
ax = axs[6]
ax.set_xlim(xgrid[1], xgrid[end])
ax.set_ylim(ygrid[1], ygrid[end])
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("attractive limit cycle"; color = "C5")

fig.subplots_adjust(bottom = 0.02, left = 0.02, top = 0.92, right = 0.97, hspace = 0.2)
wsave(plotsdir("2ddynamics"), fig)

