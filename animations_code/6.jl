using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalSystems, InteractiveChaos
using AbstractPlotting.MakieLayout

# %% Delay time impact
scene, layout = layoutscene(resolution = (1000, 800))
display(scene)
ax = layout[1, :] = LScene(scene)
sll = labelslider!(scene, "τ =", 1:100)
layout[2, :] = sll.layout
τ = sll.slider.value

ds = Systems.lorenz()
x = trajectory(ds, 100; Ttr = 100)[:, 1]

R = embed(x, 3, τ[])
Robs = Observable(R.data)

js = (206, 946)

Pobs = Observable([Point3f0(Robs[][j]) for j in js])

lines!(ax, Robs; color = COLORS[1])
scatter!(ax, Pobs; color = :red, markersize = 1000)

on(τ) do t
    Robs[] = embed(x, 3, t).data
    Pobs[] = [Point3f0(Robs[][j]) for j in js]
end
