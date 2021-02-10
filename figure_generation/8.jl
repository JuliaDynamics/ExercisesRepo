using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using InteractiveDynamics, Random
using DynamicalBilliards
import GLMakie
using GLMakie: to_color, RGBf0, RGBAf0
using DynamicalSystems

# Set style to book colors
InteractiveDynamics.obcolor(::Antidot) = to_color(COLORS[1])
InteractiveDynamics.obcolor(::Obstacle) = to_color(COLORS[1])
InteractiveDynamics.obfill(o::Antidot) = RGBAf0(0,0,0,0)
InteractiveDynamics.obls(::Antidot) = nothing


# %% Circle billiard plot
bd = Billiard(Antidot(Float32[0, 0], 1.0f0, false))

# make three particles
js = [(0.2, sind(30)), (1.5, sind(30 + sqrt(2)/10))]

colors = [to_color(COLORS[i+1]) for i in 1:length(js)]

ps = Particle{Float32}[]
for (ξ, s) in js
    x, y = from_bcoords(Float32(ξ), Float32(s), bd)
    push!(ps, Particle(x..., atan(y[2], y[1])))
end


figure, bmapax = billiard_bmap_plot(bd, ps;
colors = colors, tail = 100000, steps = 100003, backgroundcolor = RGBf0(1,1,1),
ms = 10, vr = 0.1)
GLMakie.ylims!(bmapax, 0.2, 0.9)
# GLMakie.save(joinpath(figdir, "circlebilliard.png"), figure)

# %% Chaotic billiard (sinai)
InteractiveDynamics.obcolor(::Antidot) = to_color(COLORS[1])
InteractiveDynamics.obfill(o::Antidot) = RGBAf0(0,0,0,0)
InteractiveDynamics.obls(::Antidot) = nothing

colors = [to_color(COLORS[i+1]) for i in 1:2]


bd = billiard_sinai(0.25f0, 1f0, 1f0)

ps = particlebeam(0.2f0, 0.75, π + π/4 + 0.414235, 100, 0.002)
figure, bmapax = billiard_bmap_plot(bd, ps; colors = colors,
tail = 3100, steps = 3400, backgroundcolor = RGBf0(1,1,1),
ms = 10, vr = 0.05)

# %% same code but with PyPlot for the book page
using DynamicalBilliards, PyPlot
bd = billiard_sinai() # load pre-defined billiard
ax = PyPlot.plot(bd)  # plot it
# initialize a beam of parallel particles:
ps = particlebeam(0.2, 0.75, π + π/4 + 0.414235, 100, 0.002)

for (i, p) in enumerate(ps)
    # evolve each particle for 8 collisions
    x, y = DynamicalBilliards.timeseries(p, bd, 8)
    color = (i/100, 0, 1 - i/100, 0.5)
    ax.plot(x, y; color) # plot trajectory
end



# %% Chaotic scattering
include(srcdir("style.jl"))
using PyPlot, DynamicalBilliards
fig, axs = subplots(1,3)
axs[3].axis("off")
ax3 = fig.add_subplot(1, 6, 5)
ax4 = fig.add_subplot(1, 6, 6)
DynamicalBilliards.obcolor(::Obstacle) = matplotlib.colors.to_rgb(COLORS[1])



# First plot: illustration of scattering function
# I GIVE UP, I'll Make this plot in PowerPoint
off = 0.3
r = 0.2
offset = [0.0, off]
center = 2.5
R(φ) = [cos(φ) -sin(φ);
        sin(φ)  cos(φ)]
R3 = R(2π/3)
enclosing = billiard_sinai(r, 2center, 2center)
enclosing = enclosing[2:5]

# disk0 = Disk([center,0], 1.5r, "0")
# sca(axs[1])
# plot(disk0)
# axs[1].set_xlim(-1.5 + center, 0.5 + center)
# axs[1].set_ylim(-1+0.3, 1+0.3)
axs[1].axis("off")

# p0 = Particle(center-1, 0.22, 0)
# bd0 = Billiard(enclosing..., disk0)
# x,y,vx,vy = timeseries(p0, bd0, 2; dt = 0.01)
# axs[1].plot(x[1:130],y[1:130]; color = "C1")
# x0, y0 = p0.pos
# vx0, vy0 = p0.vel
# axs[1].quiver(x0, y0, 0.08vx0, 0.08vy0; angles = "xy", scale = 0.5, width = 0.018, color="C1", zorder = 99)
# axs[1].text(x0, y0-0.2, "\$\\phi\$")
# axs[1].text(x0, y0-0.2, "\$a\$")

# Okay now clarify the scattering function
disk1 = Disk([center, center] .+ offset, r, "green")
disk2 = Disk([center, center] .+ R3*offset, r, "red")
disk3 = Disk([center, center] .+ R3*R3*offset, r, "purple")

disks = Billiard(disk1, disk2, disk3)
plot(disks; ax = axs[2])
axs[2].set_xlim(1.9,3.3)
axs[2].set_ylim(2.0,3.3)
scattering = Billiard(enclosing..., disk1, disk2, disk3)

axs[2].axis("off")

# make some particles

terminate(n, τ, i, p) = i ∈ 1:4
φs = 5π/6 .+ range(0; step = 0.001, length = 3) .+ 0.17

ps = [Particle((R(φ)*[0.5, 0] .+ [center, center])..., φ+π) for φ in φs]
for (i, p) in enumerate(ps);
    x0, y0 = p.pos
    vx0, vy0 = p.vel
    axs[2].quiver(x0, y0, 0.04vx0, 0.04vy0; angles = "xy", scale = 0.5, width = 0.01, color="C$i", zorder = 99)
    x,y = timeseries(p, scattering, terminate)
    axs[2].plot(x, y; color = "C$(i)", lw = 1, alpha = 0.75)
end


# Detailed computation of input-output function
φs = range(0, 2π/3; length = 100_000)
θs = zero(φs)

for (i, φ) in enumerate(φs);
    p = Particle((R(φ)*[0.5, 0] .+ [center, center])..., φ+π)
    timeseries!(p, scattering, terminate)
    θs[i] = atan(p.vel[2], p.vel[1])
end

ax3.plot(φs, θs .+ π; ls = "None", marker = "o", ms = 0.5, alpha = 0.5)
ax3.set_xticks([0, 2π/3])
ax3.set_xticklabels(["0", "2π/3"])
ax3.set_yticks([0, 2π])
ax3.set_yticklabels(["0", "2π"])
ax4.set_yticks([])
ax3.set_xlim(0, 2π/3)
ax4.set_ylim(0, 2π)
ax3.set_ylim(0, 2π)

φs = range(0.1555, 0.1567; length = 100_000)
θs = zero(φs)
for (i, φ) in enumerate(φs);
    p = Particle((R(φ)*[0.5, 0] .+ [center, center])..., φ+π)
    timeseries!(p, scattering, terminate)
    θs[i] = atan(p.vel[2], p.vel[1])
end
ax4.plot(φs, θs .+ π; ls = "None", marker = "o", ms = 0.5, alpha = 0.5)
ax4.set_xlim(φs[1], φs[end])
ax4.set_xticks([])
ax3.set_xlabel("\$\\phi\$", labelpad = -25)
ax3.set_ylabel("\$\\theta(\\phi)\$", labelpad = -25)
# ax3.axvspan(φs[1], φs[end]; color = "C1") # not even visible lol

fig.tight_layout(;pad = 0.25)
fig.subplots_adjust(wspace = 0.1)

wsave(plotsdir("chaoticscattering"), fig)



# %% Mushroom
w = 0.2f0
bd = billiard_mushroom(1f0, w, 1f0, 0f0)

pc = MushroomTools.randomchaotic(1f0, w, 1f0)
pc = Particle(-0.01f0, 0.2f0, sqrt(3f0)/2)

pr = Particle(0.0f0, 1.2f0, 0.0f0)
pr2 = Particle(0.0f0, 1.9f0, 0.0f0)

ps = [pc, pr, pr2]
colors = [to_color(COLORS[i+1]) for i in 1:length(ps)]

figure, bmapax = billiard_bmap_plot(bd, ps;
colors = colors, tail = 100000, steps = 1000003, backgroundcolor = RGBf0(1,1,1),
ms = 5)

# interactive_billiard_bmap(bd)

# bmapax.xticklabelsize = 32
# bmapax.yticklabelsize = 32
bmapax.xlabelsize = 40
bmapax.ylabelsize = 40
bmapax.ylabelpadding = 15
# GLMakie.save(joinpath(figdir, "mushroom.png"), figure)

# %% Natural measure of henon map
using DynamicalSystems
# using PyPlot; using3D()

ε = 0.0001
x, y, p = begin
    ds = Systems.henon()
    X = trajectory(ds, 100000000; Ttr = 100)
    p, b = ChaosTools.binhist(ε, X)
    x = [a[1] for a in b]
    y = [a[2] for a in b]
    x, y, p ./ maximum(p)
end

x ./= 2

# z = zeros(length(x))
# cmap = get_cmap("viridis")
# cs = cmap.(sqrt.(p ./ maximum(p)))
# PyPlot.figure()
# PyPlot.bar3D(x, y, z, ε/2, ε/2, p, cs)
# PyPlot.axis("off")

# GLMakie alternative:
using GLMakie
using GLMakie: Point3f0, Vec3f0, Rect3D

figure = GLMakie.meshscatter(vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(2ε, 2ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = clamp.(p, 0, 0.2))
ylims!(figure, 1.5 .* extrema(y)...)
xlims!(figure, extrema(x)...)
xticks!(figure; )
zlabel!(figure, "ρ")
display(figure)


figure
# %%

# Same but for standard map
ε = 0.01
us = [
    SVector(0.000767806, 0.00054896),
    SVector(π, 0.5),
    SVector(π, 1.0),
]
sm = Systems.standardmap(; k = 0.9)


figure = Figure()
for u in us
    x, y, p = begin
        X = trajectory(sm, 10^7, u; Ttr = 100)
        p, b = ChaosTools.binhist(ε, X)
        x = [a[1] for a in b]
        y = [a[2] for a in b]
        x, y, p ./ maximum(p)
    end
    using GLMakie

    GLMakie.meshscatter!(figure, vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(ε, ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = p, colormap = :viridis)
end

xlims!(figure, 0, 2π)
ylims!(figure, 0, 2π)
xlabel!(figure, "θ")
ylabel!(figure, "p")
zlabel!(figure, "ρ")
display(figure)


# Same but for mushroom billiard
w = 0.2
bd = billiard_mushroom(1, w, 1, 0)

pc = Particle(-0.01, 0.2, sqrt(3)/2)
pr = Particle(0.0, 1.2, 0.0)
pr2 = Particle(0.0, 1.9, 0.0)

cmaps = [:blues, :reds, :greens]
particles = [pc, pr, pr2]

ε = 0.01

figure = Figure()
for i in 1:3
    x, y, p = begin
        bmap, = boundarymap(particles[i], bd, 10^6)
        X = Dataset(bmap)
        p, b = ChaosTools.binhist(ε, X)
        x = [a[1]/2 for a in b]
        y = [a[2] for a in b]
        x, y, p ./ maximum(p)
    end

    GLMakie.meshscatter!(figure, vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(ε, ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = p, colormap = cmaps[i])
end

display(figure)

# %% Standard map island shrinking
using DynamicalSystems
grid = 0.0:0.5:2π
sm = Systems.standardmap()
ks = [0.1, 1.0, 4.0]
fig, axs = subplots(1,3; sharey = true)

for i in 1:3
    set_parameter!(sm, 1, ks[i])
    ax = axs[i]
    for θ in grid, p in grid
        u0 = SVector(θ, p)
        tr = trajectory(sm, 1000, u0)
        c = lyapunovs(sm, 4000, 1; u0 = u0)[1] > 0.01 ? "C0" : "C2"
        ax.plot(columns(tr)...; ls = "None", marker = "o", ms = 0.2, c = c)
    end
    ax.set_xlim(0, 2π)
    ax.set_ylim(0, 2π)
    ax.text(4.5, 5.5, "\$k=$(ks[i])\$"; bbox = bbox)
    ax.set_xticks(0:2:6)
end

axs[1].set_ylabel("\$p\$", labelpad = -5)
axs[2].set_xlabel("\$\\theta\$", labelpad = -15)
fig.subplots_adjust(left = 0.05, bottom = 0.14, right = 0.95, top = 0.97, wspace = 0.06)

# fsave(joinpath(figdir, "kam_sm"), fig)


# %% Recurrences
using DynamicalSystems, DynamicalBilliards
w = 0.2
pc = Particle(-0.01, 0.2, sqrt(3)/2)
pr = Particle(0.0, 1.9, 0.0)

ps = [pc, pr]
cs = ["C1", "C3"]
bd = billiard_mushroom(1, w, 1)

zbox = ((420, 130), (460, 170))

fig, axs = subplots(1, 3)
εs = [0.1, 0.1]
for i in 1:2
    ax = axs[i==1 ? 1 : 3]
    p = ps[i]
    N = 500
    bmap, = boundarymap(p, bd, N)
    tr = regularize(Dataset(bmap))
    R = RecurrenceMatrix(tr, εs[i])
    x, y = coordinates(R)
    ax.scatter(x, y, s = 1, color = cs[i])
    ax.set_aspect("equal")
    ℓ = dl_entropy(R)
    r = rt_entropy(R)
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)
    ax.set_xticks(0:250:500)
    ax.set_yticks(0:250:500)
    ax.grid()
    ax.set_title("\$H_\\ell = $(round(ℓ;digits= 1)),\\, H_r =$(round(r;digits=1))\$")
end
axs[3].set_yticklabels([])

# add zoomin
axis_zoomin!(axs[2], axs[1], zbox, zbox, "C0")
axs[2].axis("off")
axs[2].set_aspect("equal")
bmap, = boundarymap(pc, bd, 500)
tr = regularize(Dataset(bmap))
R = RecurrenceMatrix(tr, εs[1])
x, y = coordinates(R)
axs[2].scatter(x, y, s = 10, color = cs[1])
axs[2].set_xlim(zbox[1][1], zbox[2][1])
axs[2].set_ylim(zbox[1][2], zbox[2][2])

fig.subplots_adjust(left = 0.08, bottom = 0.01, right = 0.95, top = 0.97, hspace = 0.1)
# fsave(joinpath(figdir, "recurrence"), fig)
