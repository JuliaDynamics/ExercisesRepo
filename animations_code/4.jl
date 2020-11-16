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

# TODO:
# Logistic map timeseries (left panel) and f, f^2 (right panel)
# with slider for r value below them.
# %%
scene, layout = layoutscene(resolution = (1000, 800))

axts = layout[1, :] = LAxis(scene)
axmap = layout[2, :] = LAxis(scene)

rrange = 1:0.01:4.0

sll = labelslider!(scene, "r =", 0:0.01:4.0)
layout[3, :] = sll.layout
r_observable = sll.slider.value

# Timeseries plot
ns = 0:0.01:1
L = 100 # length of timeseries to be plotted
lo = Systems.logistic(0.4; r=rrange[1])
x = Observable(trajectory(lo, L))
lines!(axts, 0:L, x; color = COLORS[1], lw = 2.0)
scatter!(axts, 0:L, x; color = COLORS[1], markersize = 5)
ylims!(axts, 0, 1)

# Cobweb diagram
xs = 0:0.01:1
f1(x, r = rrange[1]) = r*x*(1-x)
f2(x, r = rrange[1]) = f1(f1(x,r), r)

f1obs = Observable(f1.(xs))
f2obs = Observable(f2.(xs))

lines!(axmap, [0,1], [0,1]; linewidth = 1, linestyle = :dash, color = COLORS[1])
f1curve = lines!(axmap, xs, f1obs; color = COLORS[2], linewidth = 2)
f2curve = lines!(axmap, xs, f2obs; color = COLORS[3], linewidth = 2)

function cobweb(t) # transform timeseries x into cobweb (point2D)
    # TODO: can be optimized to become in-place instead of allocating
    c = Point2f0[]
    for i ∈ 1:length(t)-1
        push!(c, Point2f0(t[i], t[i]))
        push!(c, Point2f0(t[i], t[i+1]))
    end
    return c
end

cobs = Observable(cobweb(x[]))
ccurve = lines!(axmap, cobs; color = COLORS[1])
cscatter = scatter!(axmap, cobs; color = COLORS[1], markersize = 2)

xlims!(axmap, 0, 1)
ylims!(axmap, 0, 1)

# On trigger r-slider update all plots:
on(r_observable) do r
    set_parameter!(lo, 1, r)
    x[] = trajectory(lo, L)
    f1obs[] = f1.(xs, r)
    f2obs[] = f2.(xs, r)
    cobs[] = cobweb(x[])
end

# Finally add three buttons to hide/show elements of the plot
f1button = LButton(scene; label = "f")
f2button = LButton(scene; label = "f²")
cbutton = LButton(scene; label = "cobweb")
layout[4, :] = buttonlayout = GridLayout(tellwidth = false)
buttonlayout[:, 1:3] = [f1button, f2button, cbutton]

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

display(scene)
