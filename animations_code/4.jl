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

figure, layout = layoutscene(resolution = (1000, 800))

ax = layout[1, :] = LAxis(figure)
sll = labelslider!(figure, "ε =", 0:0.01:1.0; sliderkw = Dict(:startvalue => 0.65))
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

display(figure)

on(ε_observable) do ε
    fs[] = dTdt.(Ts, ε)
    r, s = roots(ε)
    rootvals[] = r
    rootcols[] = s
end

# %% Logistic map cobweb and timeseries
using InteractiveChaos, Makie, DynamicalSystems

# the second range is a convenience for intermittency example of logistic
rrange = 1:0.001:4.0
# rrange = (rc = 1 + sqrt(8); [rc, rc - 1e-5, rc - 1e-3])

lo = Systems.logistic(0.4; r=rrange[1])

interactive_cobweb(
    lo, rrange, 3;
    trajcolor = COLORS[1],
    fkwargs = [(linewidth = 4.0, color = COLORS[i+1]) for i in 1:3],
)

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

# %% Cobweb for pomeaumannevile
using InteractiveChaos, Makie
using DynamicalSystems

using DynamicalSystems, PyPlot
function pm_f(x, p, n)
	γ, z = p
	s =  -2(γ+1)
    if x < -0.5
		s*(x+0.5) - γ
        # -4x - 3
    elseif -0.5 ≤ x ≤ 0.5
        @inbounds γ*x*(1 + abs(2x)^(z-1))
    else
        s*(x-0.5) + γ
    end
end

rs = 0.7:0.01:1
ds = DiscreteDynamicalSystem(pm_f, 0.3, [rs[1], 2.6])

interactive_cobweb(
    ds, rs, 1;
    trajcolor = COLORS[1],
    fkwargs = [(linewidth = 4.0, color = COLORS[i+1]) for i in 1:3],
	xmin = -1,
	xmax = 1
)
