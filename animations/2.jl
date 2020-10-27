using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))

# %% Animation of changing ε in curve of 1D energy balance model
using Makie, DynamicalSystems, InteractiveChaos
using AbstractPlotting.MakieLayout

function fitzhugh(x, p, t)
    u, w = x
    I, a, b = p
    udot = a*u*(u-b)*(1-u) - w + I
    wdot = 0.01*(u-w)
    return SVector(udot, wdot)
end

a = 3.0
b = 0.2
I = 0.0
p = [I, a, b]
fh = ContinuousDynamicalSystem(fitzhugh, rand(2), p)

# Layout plot:
scene, layout = layoutscene(resolution = (1000, 800))
ax = layout[1,1] = LAxis(scene)
display(scene)
ax.xlabel = "u"
ax.ylabel = "w"

# plot the two nullclines
us = -0.4:0.01:1.2
w1 = us
w2 = @. a*us*(us-b)*(1-us) + I
lines!(ax, us, w1; color = COLORS[1], linewidth = 2, linestyle = :dash)
lines!(ax, us, w2; color = COLORS[2], linewidth = 2, linestyle = :dot)

spoint = select_point(ax.scene)
on(spoint) do pos
    tr = trajectory(fh, 10000, SVector(pos...))
    lines!(ax.scene, columns(tr)...; color = InteractiveChaos.randomcolor())
end
