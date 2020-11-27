using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using Makie, DynamicalSystems, InteractiveChaos
using AbstractPlotting.MakieLayout
using Roots

# %% Interactive bifurcations for 1D energy balance

αtan(T) = 0.5 - (0.4/π)*atan((T-263)/2)
dTdt(T, ε = 0.65, α=αtan, s= 1.0) = s*(1 - α(T)) - 1.6e-10 * ε * T^4

Ts = 200:0.5:320.0

function roots(ε)
    f = T -> dTdt(T, ε)
    r = Roots.find_zeros(f, Ts[1], Ts[end])
    s = [to_color(COLORS[isodd(i) ? 4 : 3]) for i in 1:length(r)]
    return [Point2f0(x, 0) for x in r], s # roots, stability
end

scene, layout = layoutscene(resolution = (1000, 800))

ax = layout[1, :] = LAxis(scene)
sll = labelslider!(scene, "ε =", 0:0.01:1.0; sliderkw = Dict(:startvalue => 0.65))
layout[2, :] = sll.layout
ε_observable = sll.slider.value

# initialize plot
fs = Observable(dTdt.(Ts))
lines!(ax, Ts, fs; color = COLORS[2], linewidth = 4)
lines!(ax, Ts, zero(Ts); color = COLORS[1])
r, s = roots(0.65)
rootvals = Observable(r)
rootcols = Observable(s)
scatter!(ax, rootvals; color = rootcols, markersize = 10)

display(scene)

on(ε_observable) do ε
    fs[] = dTdt.(Ts, ε)
    r, s = roots(ε)
    rootvals[] = r
    rootcols[] = s
end

# %% Logistic map cobweb and timeseries
scene, layout = layoutscene(resolution = (1000, 800))

axts = layout[1, :] = LAxis(scene)
axmap = layout[2, :] = LAxis(scene)

# `rrange` decides the r-slider values.
# the second range is a convenience for intermittency
rrange = 1:0.001:4.0
# rrange = (rc = 1 + sqrt(8); [rc, rc - 1e-5, rc - 1e-2])

slr = labelslider!(scene, "r =", rrange)
layout[3, :] = slr.layout
r_observable = slr.slider.value

sln = labelslider!(scene, "n =", 10:1000)
layout[4, :] = sln.layout
L = sln.slider.value

# Timeseries plot
function seriespoints(x)
    n = 0:length(x)+1
    c = [Point2f0(n[i], x[i]) for i in 1:length(x)]
end

lo = Systems.logistic(0.4; r=rrange[1])
x = Observable(trajectory(lo, L[]))
xn = lift(a -> seriespoints(a), x)
lines!(axts, xn; color = COLORS[1], lw = 2.0)
scatter!(axts, xn; color = COLORS[1], markersize = 5)
ylims!(axts, 0, 1)

# Cobweb diagram
xs = 0:0.001:1
f1(x, r = rrange[1]) = r*x*(1-x)
f2(x, r = rrange[1]) = f1(f1(x,r), r)
f3(x, r = rrange[1]) = f2(f1(x,r), r)

f1obs = Observable(f1.(xs))
f2obs = Observable(f2.(xs))
f3obs = Observable(f3.(xs))

lines!(axmap, [0,1], [0,1]; linewidth = 1, linestyle = :dash, color = COLORS[1])
f1curve = lines!(axmap, xs, f1obs; color = COLORS[2], linewidth = 2)
f2curve = lines!(axmap, xs, f2obs; color = COLORS[3], linewidth = 2)
f3curve = lines!(axmap, xs, f3obs; color = COLORS[4], linewidth = 2)

function cobweb(t) # transform timeseries x into cobweb (point2D)
    # TODO: can be optimized to become in-place instead of allocating
    c = Point2f0[]
    for i ∈ 1:length(t)-1
        push!(c, Point2f0(t[i], t[i]))
        push!(c, Point2f0(t[i], t[i+1]))
    end
    return c
end

cobs = lift(a -> cobweb(a), x)
ccurve = lines!(axmap, cobs; color = COLORS[1])
cscatter = scatter!(axmap, cobs; color = COLORS[1], markersize = 2)

xlims!(axmap, 0, 1)
ylims!(axmap, 0, 1)

# On trigger r-slider update all plots:
on(r_observable) do r
    set_parameter!(lo, 1, r)
    x[] = trajectory(lo, L[])
end

on(L) do l
    x[] = trajectory(lo, l)
    xlims!(axts, 0, l)
end

# Finally add buttons to hide/show elements of the plot
f1button = LButton(scene; label = "f")
f2button = LButton(scene; label = "f²")
f3button = LButton(scene; label = "f³")
cbutton = LButton(scene; label = "cobweb")
layout[5, :] = buttonlayout = GridLayout(tellwidth = false)
buttonlayout[:, 1:4] = [f1button, f2button, f3button, cbutton]

# And add triggering for buttons
on(cbutton.clicks) do click
    ccurve.attributes.visible = !(ccurve.attributes.visible[])
    cscatter.attributes.visible = !(cscatter.attributes.visible[])
end
on(f1button.clicks) do click
    f1curve.attributes.visible = !(f1curve.attributes.visible[])
end
on(f2button.clicks) do click
    f2curve.attributes.visible = !(f2curve.attributes.visible[])
end
on(f3button.clicks) do click
    f3curve.attributes.visible = !(f3curve.attributes.visible[])
end

display(scene)

# %% Interactive orbit diagram for logistic map
using InteractiveChaos, Makie
using DynamicalSystems

ds = Systems.logistic()
p_min, p_max = 1.0, 4.0
t = "orbit diagram for the logistic map"

oddata = interactive_orbitdiagram(
    ds, 1, p_min, p_max, 1;
    parname = "r", title = t
)
