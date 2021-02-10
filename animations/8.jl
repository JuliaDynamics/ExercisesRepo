using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using GLMakie, DynamicalBilliards, InteractiveChaos

# Set style to book colors
InteractiveChaos.obcolor(::Antidot) = to_color(COLORS[1])
InteractiveChaos.obcolor(::Obstacle) = to_color(COLORS[1])
InteractiveChaos.obfill(o::Antidot) = RGBAf0(0,0,0,0)
InteractiveChaos.obls(::Antidot) = nothing

# %% Circle billiard animation
bd = Billiard(Antidot([0.0, 0.0], 1.0, false))

js = [(0.2, sind(30)), (1.5, sind(30 + sqrt(2)/10))]
colors = [to_color(COLORS[i+1]) for i in 1:length(js)]

ps = Particle{Float64}[]
for (ξ, s) in js
    x, y = from_bcoords(ξ, s, bd)
    push!(ps, Particle(x..., atan(y[2], y[1])))
end

interactive_billiard(
    bd, ps;
    colors, tail = 20000, backgroundcolor = RGBf0(1,1,1),
)

# %% Present it as a boundary map plot
billiard_bmap_plot(
    bd, ps;
    colors, tail = 20000, backgroundcolor = RGBf0(1,1,1),
    steps = 10000, dt = 0.01,
)

# %% Sinai billiard
bd = billiard_sinai(0.25)
ps = particlebeam(0.2, 0.75, π + π/4 + 0.414235, 100, 0.002)

interactive_billiard(
    bd, ps;
    tail = 500, backgroundcolor = RGBf0(1,1,1),
    vr = 0.05
)

# %% Mixed state space example
w = 0.2f0
bd = billiard_mushroom(1f0, w, 1f0, 0f0)

pc = MushroomTools.randomchaotic(1f0, w, 1f0)
pc = Particle(-0.01f0, 0.2f0, sqrt(3f0)/2)

pr = Particle(0.0f0, 1.2f0, 0.0f0)
pr2 = Particle(0.0f0, 1.9f0, 0.0f0)

ps = [pc, pr, pr2]
colors = [to_color(COLORS[i+1]) for i in 1:length(ps)]

figure, bmapax = interactive_billiard(bd, ps;
colors = colors, tail = 10000, )

# %% Zooming into standard map
using DynamicalSystems
figure = Figure(resolution = (1000,1000))
ax = figure[1,1] = Axis(figure)
deactivate_interaction!(ax, :rectanglezoom)
rect = select_rectangle(ax.scene)
sm = Systems.standardmap(; k = 1.0)
g = 10 # grid density
inrect(u, r) = (r.origin[1] ≤ u[1] ≤ r.origin[1]+r.widths[1]) && (r.origin[2] ≤ u[2] ≤ r.origin[2] + r.widths[2])
integ = integrator(sm)
marker = InteractiveChaos.MARKER

N = 100 # how many points to plot inside the rect
Nslider = labelslider!(
    figure, "N =", round.(Int, range(100, 100000; length = 100))
)
N = Nslider.slider.value

recompute = Button(figure; label = "recompute")

figure[2,1][1,1] = recompute
figure[2,1][1,2] = Nslider.layout

sco = Observable([Point2f0(sm.u0)])

scatter!(ax, sco; markersize = 1*px, markerstrokewidth = 0)

on(rect) do r # new subgrid selected
    xs = range(r.origin[1], r.origin[1]+r.widths[1]; length = g)
    ys = range(r.origin[2], r.origin[2]+r.widths[2]; length = g)
    s = sco[]
    n = N[]
    empty!(s)
    for x in xs
        for y in ys
            DynamicalSystems.reinit!(integ, SVector(x, y))
            push!(s, integ.u)
            i = 0
            for i in 1:n
                DynamicalSystems.step!(integ)
                if inrect(integ.u, r)
                    push!(s, integ.u)
                end
            end
        end
    end
    sco[] = s
    xlims!(ax, xs[1], xs[end])
    ylims!(ax, ys[1], ys[end])
end

on(recompute.clicks) do c
    rect[] = rect[]
end

rect[] = FRect2D([0,0],[2π,2π])
ax.xticklabelsvisible = false
ax.yticklabelsvisible = false

display(figure)

# %% Ergodicity in the Sinai billiard
bd = billiard_sinai(0.25)

interactive_billiard_bmap(bd)
