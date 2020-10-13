# Animations and extra plots used in chapter 1
using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))

# Animate evolution of trajectories in Lorenz
using InteractiveChaos
using DynamicalSystems, Makie
using OrdinaryDiffEq

ds = Systems.lorenz()

u1 = [10,20,40.0]
u2 = [10,20,40.0 + 1e-3]
u3 = [20,10,40.0]
u0s = [u1, u2, u3]

diffeq = (alg = Tsit5(), dtmax = 0.01)

scene, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 10000, diffeq, colors = COLORSCHEME,
)

# %% Add poincare plane
using AbstractPlotting.MakieLayout

o = Point3f0(-25, 0, 0)
w = Point3f0(50, 0, 50)
p = FRect3D(o, w)
# mesh!(main, p; color = RGBAf0(0.2, 0, 0.8), transparent = true)

# These lines are necessary for transparent planes
a = RGBAf0(0,0,0,0)
c = RGBAf0(0.2, 0.2, 0.8, 1.0)
img = AbstractPlotting.ImagePattern([c a; a c]); # This throws an error if it shows, you can ignore that
mesh!(main, p; color = img)

# Plot Poincare sos
psosplot = layout[:, 2] = LAxis(scene)
psos = poincaresos(ds, (2, 0.0), 2000.0)
scatter!(psosplot, psos[:, 1], psos[:, 3])

display(scene)

# %%
ds = Systems.henonheiles()

u0s = [[0.0, -0.25, 0.42081, 0.0],
[0.0, 0.1, 0.5, 0.0],
[0.0, -0.31596, 0.354461, 0.0591255]]

diffeq = (alg = Vern9(), dtmax = 0.01)
idxs = (1, 2, 4)
colors = ["#233B43", "#499cbf", "#E84646"]

scene, main, layout, obs = interactive_evolution(
    ds, u0s; idxs, tail = 2000, diffeq, colors
)
main.scene[Axis][:names, :axisnames] = ("q₁", "q₂", "p₂")
