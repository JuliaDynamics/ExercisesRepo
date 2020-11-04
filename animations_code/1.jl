using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))

# %% The very first code snippet of the course
using DynamicalSystems # load the library

function lorenz_rule(u, p, t)
    σ, ρ, β = p
    x, y, z = u
    dx = σ*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - β*z
    return SVector(dx, dy, dz) # Static Vector
end

p  = [10.0, 28.0, 8/3] # parameters: σ, ρ, β
u₀ = [0, 10.0, 0]      # initial state
# create an instance of a `DynamicalSystem`
lorenz = ContinuousDynamicalSystem(lorenz_rule, u₀, p)

T  = 100.0 # total time
dt = 0.01  # sampling time
A  = trajectory(lorenz, T; dt)

# %% Animate evolution of trajectories in the Standard map
using InteractiveChaos
using DynamicalSystems
import Makie
using OrdinaryDiffEq

# Standard map trajectories
ds = Systems.standardmap(; k = 1.0)
u0s = [[θ, p] for θ ∈ 0:2π for p ∈ 0:2π]
lims = ((0, 2π), (0, 2π))

scene, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 1000, lims,
)
main.xlabel = "θ"
main.ylabel = "p"

# %% Lorenz system trajectories
ds = Systems.lorenz()

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
