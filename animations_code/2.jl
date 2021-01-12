using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using Makie, DynamicalSystems, InteractiveChaos
using AbstractPlotting.MakieLayout

# %% Fitzhugh clickable phase space
a = 3.
b = -0.05
ε = 0.01
I = 0.0

fh = Systems.fitzhugh_nagumo(zeros(2); a, b, ε, I)

# Layout plot:
figure, layout = layoutscene(resolution = (1000, 800))
ax = layout[1,1] = LAxis(figure)
display(figure)
ax.xlabel = "u"
ax.ylabel = "w"
ax.title = "FitzHugh Nagumo, a=$a, b=$b"

# plot the two nullclines
us = -0.4:0.01:1.2
w1 = us
w2 = @. a*us*(us-b)*(1-us) + I
lines!(ax, us, w1; color = COLORS[1], linewidth = 2, linestyle = :dash)
lines!(ax, us, w2; color = COLORS[2], linewidth = 2, linestyle = :dot)

# Create trajectories on point selection
spoint = select_point(ax.figure)
on(spoint) do pos
    tr = trajectory(fh, 10000, SVector(pos...))
    lines!(ax, columns(tr)...;
        color = InteractiveChaos.randomcolor(), linewidth = 4.0
    )
end

# %% orbit evolution in a 2D torus
using OrdinaryDiffEq

R = 2.0
r = 1.0

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
u0s = [[0, 0, ω] for ω ∈ [sqrt(7), 3]]
diffeq = (alg = Tsit5(), dtmax = 0.01)

lims = ((-R-r, R+r), (-R-r, R+r), (-3r, 3r))

plotkwargs = [(linewidth = i^2*1.0,) for i in 1:length(u0s)]

figure, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 20000, diffeq, colors = COLORS, transform = torus, lims,
    plotkwargs
)
