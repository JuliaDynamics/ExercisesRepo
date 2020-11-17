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


# %% Magnetic pendulum
using LinearAlgebra
ma = Systems.magnetic_pendulum(α=0.2, ω=1.0, d=0.3)
∫ = integrator(ma)

g = range(-5, 5; length = 1500)
c = fill(0, length(g), length(g))
t = similar(c, Float64)
@show length(c)

for (i, x) ∈ enumerate(g)
    @show x
    for (j, y) in enumerate(g)
    reinit!(∫, SVector(x, y, 0, 0))
    t0 = ∫.t0
    step!(∫, 100.0)
    while ∫.u[3]^2 + ∫.u[4]^2 > 1e-3
        step!(∫)
    end
    s = SVector(∫.u[1], ∫.u[2])
    k = findmin([(s-m)⋅(s-m) for m in ma.f.magnets])[2]
    # scatter(x, y, color = "C$(k-1)", s = 1)
    t[i,j] = ∫.t - ∫.t0
    c[i,j] = k
    # break
end
end

cd(@__DIR__)
using FileIO
save("magneticpendulum.bson", Dict(:c => c, :g => g, :t => t))

# %% Magnetic pendulum zoomed
using LinearAlgebra
ma = Systems.magnetic_pendulum(α=0.2, ω=1.0, d=0.3)
∫ = integrator(ma)

gx = range(1.80, 1.95; length = 1000)
gy = range(0, 0.12; length = 1000)
c = fill(0, length(gx), length(gy))
t = similar(c, Float64)
@show length(c)

for (i, x) ∈ enumerate(gx)
    @show x
    for (j, y) in enumerate(gy)
        reinit!(∫, SVector(x, y, 0, 0))
        t0 = ∫.t0
        step!(∫, 100.0)
        while ∫.u[3]^2 + ∫.u[4]^2 > 1e-3
            step!(∫)
        end
        s = SVector(∫.u[1], ∫.u[2])
        k = findmin([(s-m)⋅(s-m) for m in ma.f.magnets])[2]
        # scatter(x, y, color = "C$(k-1)", s = 1)
        t[i,j] = ∫.t - ∫.t0
        c[i,j] = k
        # break
    end
end

cd(@__DIR__)
using FileIO
save("magneticpendulum_zoom.bson", Dict(:c => c, :gx => gx, :gy => gy, :t => t))


# %% Plot this
ma = Systems.magnetic_pendulum(α=0.2, ω=1.0, d=0.3)

cd(@__DIR__)
using FileIO
a = load("magneticpendulum.bson")
g = a[:g]; c = a[:c]

LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])

fig = figure(figsize=(figx/2, figx/2))
pcolormesh(g, g, c'; cmap = cmap, shading = "gouraud")
gca().set_aspect("equal")
xticks([-5, 5])
yticks([-5, 5])
xlabel("\$x\$", labelpad=-30)
ylabel("\$y\$", labelpad=-30)
for m in ma.f.magnets
    scatter(m[1], m[2]; color = "white", edgecolor = "black", zorder = 99, s = 100)
end
tight_layout()
subplots_adjust(bottom = 0.1, top = 0.95)


# plot rando initial condition
# tr = trajectory(ma, 10000, [-1.91228, -3.59712, 0, 0]; dt = 0.1)
# x, y = columns(tr)
# plot(x, y, "C3", lw = 1, alpha = 0.5, zorder = 1)
fsave(joinpath(figdir, "magneticpendulum"), fig)


az = load("magneticpendulum_zoom.bson")
gx = az[:gx]; gy = az[:gy]; c = az[:c]

LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])

fig = figure(figsize=(figx/2, figx/2))
pcolormesh(gx, gy, c'; cmap = cmap, shading = "gouraud")
gca().set_aspect("equal")

xticks([1.8, 1.95], size = 20)
yticks([0, 0.12], size = 20)

fig.savefig(joinpath(figdir, "magneticpendulum_zoom.png"))

# xticks([-5, 5])
# yticks([-5, 5])
