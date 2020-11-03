using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using Makie, DynamicalSystems, InteractiveChaos
using AbstractPlotting.MakieLayout

# %% Animation of changing ε in curve of 1D energy balance model
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

# %% orbit evolution in a 2D torus
using OrdinaryDiffEq

R = 2.0
r = 0.5

function torus(u)
    θ, φ = u
    x = (R + r*cos(θ))*cos(φ)
    y = (R + r*cos(θ))*sin(φ)
    z = r*sin(θ)
    return SVector(x, y, z)
end

function quasiperiodic_f(u, p, t)
    # here we make the frequency ratio a state variable, because the
    # trajectory_animator application doesn't allow different parameters for different
    # initial conditions
    θ, φ, ω = u
    θdot = ω
    φdot = 1.0
    return SVector(θdot, φdot, 0.0)
end

ds = ContinuousDynamicalSystem(quasiperiodic_f, [0.0, 0.0, 0.0], nothing)
u0s = [[0, 0, ω] for ω ∈ [3, sqrt(7), 5]]
diffeq = (alg = Tsit5(), dtmax = 0.01)

lims = ((-R-r, R+r), (-R-r, R+r), (-3r, 3r))

scene, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 20000, diffeq, colors = COLORS, transform = torus, lims,
    plotkwargs = (linewidth = 1.0,),
)
