include("style.jl")
using InteractiveChaos, Random
using DynamicalBilliards
import Makie
using Makie: to_color, RGBf0, RGBAf0
using DynamicalSystems

# Set style to book colors
InteractiveChaos.obcolor(::Antidot) = to_color(COLORS[1])
InteractiveChaos.obcolor(::Obstacle) = to_color(COLORS[1])
InteractiveChaos.obfill(o::Antidot) = RGBAf0(0,0,0,0)
InteractiveChaos.obls(::Antidot) = nothing


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


scene, bmapax = billiard_bmap_plot(bd, ps;
colors = colors, tail = 100000, steps = 100003, backgroundcolor = RGBf0(1,1,1),
ms = 10, vr = 0.1)
Makie.ylims!(bmapax, 0.2, 0.9)
# Makie.save(joinpath(figdir, "circlebilliard.png"), scene)

# %% Chaotic billiard (sinai)
InteractiveChaos.obcolor(::Antidot) = to_color(COLORS[1])
InteractiveChaos.obfill(o::Antidot) = RGBAf0(0,0,0,0)
InteractiveChaos.obls(::Antidot) = nothing

colors = [to_color(COLORS[i+1]) for i in 1:2]


bd = billiard_sinai(0.25f0, 1f0, 1f0)

ps = particlebeam(0.2f0, 0.75, π + π/4 + 0.414235, 100, 0.002)
scene, bmapax = billiard_bmap_plot(bd, ps; colors = colors,
tail = 3100, steps = 3400, backgroundcolor = RGBf0(1,1,1),
ms = 10, vr = 0.05)

# %% Mushroom
w = 0.2f0
bd = billiard_mushroom(1f0, w, 1f0, 0f0)

pc = MushroomTools.randomchaotic(1f0, w, 1f0)
pc = Particle(-0.01f0, 0.2f0, sqrt(3f0)/2)

pr = Particle(0.0f0, 1.2f0, 0.0f0)
pr2 = Particle(0.0f0, 1.9f0, 0.0f0)

ps = [pc, pr, pr2]
colors = [to_color(COLORS[i+1]) for i in 1:length(ps)]

scene, bmapax = billiard_bmap_plot(bd, ps;
colors = colors, tail = 100000, steps = 1000003, backgroundcolor = RGBf0(1,1,1),
ms = 5)

# interactive_billiard_bmap(bd)

# bmapax.xticklabelsize = 32
# bmapax.yticklabelsize = 32
bmapax.xlabelsize = 40
bmapax.ylabelsize = 40
bmapax.ylabelpadding = 15
# Makie.save(joinpath(figdir, "mushroom.png"), scene)

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

# Makie alternative:
using Makie
using Makie: Point3f0, Vec3f0, Rect3D

scene = Makie.meshscatter(vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(2ε, 2ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = clamp.(p, 0, 0.2))
ylims!(scene, 1.5 .* extrema(y)...)
xlims!(scene, extrema(x)...)
xticks!(scene; )
zlabel!(scene, "ρ")
display(scene)


scene
# %%

# Same but for standard map
ε = 0.01
us = [
    SVector(0.000767806, 0.00054896),
    SVector(π, 0.5),
    SVector(π, 1.0),
]
sm = Systems.standardmap(; k = 0.9)


scene = Scene()
for u in us
    x, y, p = begin
        X = trajectory(sm, 10^7, u; Ttr = 100)
        p, b = ChaosTools.binhist(ε, X)
        x = [a[1] for a in b]
        y = [a[2] for a in b]
        x, y, p ./ maximum(p)
    end
    using Makie

    Makie.meshscatter!(scene, vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(ε, ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = p, colormap = :viridis)
end

xlims!(scene, 0, 2π)
ylims!(scene, 0, 2π)
xlabel!(scene, "θ")
ylabel!(scene, "p")
zlabel!(scene, "ρ")
display(scene)


# Same but for mushroom billiard
w = 0.2
bd = billiard_mushroom(1, w, 1, 0)

pc = Particle(-0.01, 0.2, sqrt(3)/2)
pr = Particle(0.0, 1.2, 0.0)
pr2 = Particle(0.0, 1.9, 0.0)

cmaps = [:blues, :reds, :greens]
particles = [pc, pr, pr2]

ε = 0.01

scene = Scene()
for i in 1:3
    x, y, p = begin
        bmap, = boundarymap(particles[i], bd, 10^6)
        X = Dataset(bmap)
        p, b = ChaosTools.binhist(ε, X)
        x = [a[1]/2 for a in b]
        y = [a[2] for a in b]
        x, y, p ./ maximum(p)
    end

    Makie.meshscatter!(scene, vec(Point3f0.(x, y, 0.0)),
    markersize=Vec3f0.(ε, ε, p), marker=Rect3D(Vec3f0(0), Vec3f0(1)),
    limits=Rect3D(Vec3f0(0), Vec3f0(1)),
    color = p, colormap = cmaps[i])
end

display(scene)

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
