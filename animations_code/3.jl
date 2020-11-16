using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using Makie, DynamicalSystems, InteractiveChaos
using AbstractPlotting.MakieLayout

using OrdinaryDiffEq: Tsit5, Vern9
ds = Systems.lorenz()
u0 = 3ones(dimension(ds))
u0s =  [u0 .+ 1e-3i for i in 1:3]

diffeq = (alg = Vern9(), dtmax = 0.01)

scene, main, layout, obs = interactive_evolution(
    ds, u0s; tail = 100, diffeq, colors = COLORS,
)


# %% Interactive GALI psos for henon heiles
# TODO: I need to fix PSOS to not autozoom.

using InteractiveChaos, Makie, OrdinaryDiffEq, DynamicalSystems
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

hh = Systems.henonheiles()

potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
E = energy(get_state(hh)...)

function complete(y, py, x)
    V = potential(x, y)
    Ky = 0.5*(py^2)
    Ky + V ≥ E && error("Point has more energy!")
    px = sqrt(2(E - V - Ky))
    ic = [x, y, px, py]
    return ic
end

plane = (1, 0.0) # first variable crossing 0

cmap = cgrad(:viridis)
function galicolor(u)
    g, t = gali(hh, 4000, 4; u0 = u)
    @show t[end]
    v = clamp(t[end]/500, 0, 1)
    # return cmap.colors[round(Int, v*255 + 1)]
    RGBf0(0, 0, clamp(v, 0, 1))
end

state, scene = interactive_poincaresos(
    hh, plane, (2, 4), complete;
    labels = ("q₂" , "p₂"), color = galicolor, diffeq...
)

# %% 0-1 test animation
