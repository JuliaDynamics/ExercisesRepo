# Animations and extra plots used in chapter 1
using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))

# Animate evolution of trajectories in Lorenz
using InteractiveChaos
using DynamicalSystems
using Makie
using OrdinaryDiffEq

# Standard map trajectories
ds = Systems.standardmap()
u0s = [[θ, p] for θ ∈ 0:2π for p ∈ 0:2π]
lims = ((0, 2π), (0, 2π))

scene, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 100000, lims,
)
main.xlabel = "θ"
main.ylabel = "p"

# %% Lorenz system trajectories
ds = Systems.lorenz()

u1 = [10,20,40.0]
u2 = [10,20,40.0 + 1e-3]
u3 = [20,10,40.0]
u0s = [u1, u2, u3]
u0s =  [[10,20,40.0] .+ rand(3) for _ in 1:7]

diffeq = (alg = Tsit5(), dtmax = 0.01)

scene, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 1000, diffeq, colors = COLORS,
)

# %% Add poincare plane
using AbstractPlotting.MakieLayout

o = Point3f0(-25, 0, 0)
w = Point3f0(50, 0, 50)
p = FRect3D(o, w)

# These lines are necessary for transparent planes
a = RGBAf0(0,0,0,0)
c = RGBAf0(0.2, 0.2, 0.8, 1.0)
img = AbstractPlotting.ImagePattern([c a; a c]);
mesh!(main, p; color = img)

# %% Plot Poincare sos
psosplot = layout[:, 2] = LAxis(scene)
psos = poincaresos(ds, (2, 0.0), 2000.0)
scatter!(psosplot, psos[:, 1], psos[:, 3])

display(scene)

# %% Henon heiles system trajectories
ds = Systems.henonheiles()

u0s = [[0.0, -0.25, 0.42081, 0.0],
[0.0, 0.1, 0.5, 0.0],
[0.0, -0.31596, 0.354461, 0.0591255]]

diffeq = (alg = Vern9(), dtmax = 0.01)
idxs = (1, 2, 4)

scene, main, layout, obs = interactive_evolution(
    ds, u0s; idxs, tail = 2000, diffeq, colors = COLORS,
)
main.scene[Axis][:names, :axisnames] = ("q₁", "q₂", "p₂")
