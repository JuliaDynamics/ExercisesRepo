# Logistic map orbit diagram
using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalSystems, InteractiveDynamics

ds = Systems.logistic()
p_min, p_max = 1.0, 4.0
t = "orbit diagram for the logistic map"

fig, oddata = interactive_orbitdiagram(
    ds, 1, p_min, p_max, 1;
    parname = "r", title = t
)

record_interaction(
    fig, projectdir("animations", "4", "logistic_orbitdiagram.mp4"); 
    total_time=15, sleep_time = 2
)
