using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))

using InteractiveChaos
using DynamicalSystems
using GLMakie

# Standard map trajectories
ds = Systems.standardmap(; k = 1.0)
u0s = [[θ, p] for θ ∈ 0:2π for p ∈ 0:2π]
lims = ((0, 2π), (0, 2π))

figure, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 1000, lims,
)
main.xlabel = "θ"
main.ylabel = "p"
