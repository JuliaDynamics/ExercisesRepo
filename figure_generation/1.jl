using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# %% phase space plots
fig = figure(figsize = (figx, figx/2))

ax = subplot(231)

ds = Systems.standardmap(k=1)
ic = [[0.1, 0.1], [2.5, 0.4], [1.88, 3.25]]
xs = range(0, stop = 2π, length = 7);
ys = range(0, stop = 2π, length = 3);
ys = copy(xs)
iters = 2000
dataset = Dataset{2, Float64}()
for (i, c) in enumerate(ic)
    tr = trajectory(ds, iters, c)
    ax.scatter(columns(tr)..., s = 5)
end
# ax.set_xticks([0, 2π])
ax.set_xticklabels([])
ax.set_yticklabels([])
# ax.set_yticks([0, 2π])
# ax.set_yticks([0, 2π])
ax.set_yticklabels([])
ax.set_xlabel("\$\\theta\$", labelpad = -5)
ax.set_xlim(0, 2π)
ax.set_ylim(0, 2π)
ax.set_ylabel("\$p\$",labelpad=0)
ax.set_title("Standard map")

ax = subplot(234)
n = 100
tr = trajectory(ds, n, ic[1])
ax.plot(0:n, tr[:, 2], c = COLORS[1], marker = "o", lw = 1.0)
ax.set_yticks([0, 2π])
ax.set_yticklabels(["0", "2\$\\pi\$"])
ax.set_ylabel("\$p\$",labelpad=-15)
ax.set_xticks([0, n])
ax.set_xlabel("\$n\$", labelpad = -15)

ds = Systems.lorenz()
using3D()
ax = subplot(232, projection = "3d")

ax.set_title("Lorenz-63", pad = 24)

llw = 2.0
tr = trajectory(ds, 100; Ttr = 100, dt = 0.01)
ax.plot3D(columns(tr)..., color = COLORS[1], lw = llw)
tr = trajectory(ds, 10; Ttr = 1000, dt = 0.01)
ax.plot3D(columns(tr)..., color = COLORS[2], lw = llw)
tr = trajectory(ds, 10, [10,20,40.0]; Ttr = 10, dt = 0.01)
ax.plot3D(columns(tr)..., color = COLORS[3], lw = llw)
ax.set_xlabel("\$x\$", labelpad = 0)
ax.set_ylabel("\$y\$", labelpad = 0)
ax.set_zlabel("\$z\$", labelpad = 0)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
# ax.set_xticks([-20,20])
# ax.set_yticks([-20,20])
# ax.set_zticks([10,40])

tn = 30
ax = subplot(235)
tr = trajectory(ds, tn; Ttr = 40, dt = 0.01)
ax.plot(0:0.01:tn, tr[:, 1], c = COLORS[1])

ax.set_yticks([-15,15])
ax.set_ylabel("\$x\$", labelpad = -30)
ax.set_xticks([0, tn])
ax.set_xlabel("\$t\$", labelpad = -15)

ds = Systems.henonheiles()

u0s = (
    [0.0, -0.25, 0.42, 0.0], # chaotic
    [0.0, 0.1, 0.5, 0.0], # quasiperiodic
    [0.0, 0.30266571044921875, 0.4205654433900762, 0.0], # periodic
)

ax = subplot(233, projection = "3d")
ax.clear()
ax.set_title("Hénon–Heiles", pad = 24)
for (i, u0) in enumerate(u0s)
    tr = trajectory(ds, 50, u0; Ttr = 50)
    ax.plot3D(tr[:, 1], tr[:, 2], tr[:, 3], c = COLORS[i], alpha = 1)
end
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_xlabel("\$x\$")
ax.set_ylabel("\$y\$")
ax.set_zlabel("\$p_x\$",labelpad=-5)

ax = subplot(236)
n = 150
tr = trajectory(ds, n, [0.0, -0.25, 0.48, 0.0]; Ttr = 100)
ax.plot(0:0.01:n, tr[:, 2], c = COLORS[1])
ax.set_xticks([0, n])
ax.set_yticks([-0.5, 0.5])
ax.set_ylabel("\$y\$", labelpad = -30)
ax.set_xlabel("\$t\$", labelpad = -15)


fig.tight_layout()
fig.subplots_adjust(top = 0.93, wspace = 0.4, hspace = 0.25, bottom = 0.1, left = 0.05)
wsave(plotsdir("trajectories"), fig)
