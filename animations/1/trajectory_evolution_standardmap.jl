using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))

using InteractiveDynamics
using DynamicalSystems
using GLMakie

# Standard map trajectories
ds = Systems.standardmap(; k = 1.0)
u0s = [[θ, p] for θ ∈ 0:2π for p ∈ 0:2π]
lims = ((0, 2π), (0, 2π))

figure, obs = interactive_evolution(
    ds, u0s; tail = 1000, lims,
)

# extract the main axis to add custom labels:
main = content(figure[1, 1])
main.xlabel = "θ"
main.ylabel = "p"
