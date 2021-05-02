using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Statistics, DelimitedFiles

# %% Carbon cycle extinction events
# model of Rothman, 10.1073/pnas.1905164116
# main dynamic rule, eqs.(11, 12) of paper
function excitable_carbon(u, p, t)
    c, w = u
    μ, b, θ, cₓ, cₚ, ν, w₀, γ, f₀, c_f, β = p
    s = excitable_carbon_s(c, cₚ, γ)
    sbar = 1 - excitable_carbon_s(c, cₓ, γ)
    f = excitable_carbon_f(c, f₀, β, c_f)
    cdot = μ*(1 - b*s - θ*sbar - ν) + w - w₀
    cdot *= f
    wdot = μ*(1 - b*s + θ*sbar + ν) - w + w₀
    return SVector(cdot, wdot)
end
# functions used in equations of motion
excitable_carbon_f(c, f₀, β, c_f) = f₀*(c^β / (c^β + c_f^β))
excitable_carbon_s(c, cₚ, γ) = c^γ/(c^γ + cₚ^γ)

u0 = [100.0, 2800.0] # random, from Fig. 4

# Parameter default values in Table 1 of appendix, but cₓ same as Fig. 5
p0 = [250, 4, 5, 55, 110, 0, 2000, 4, 0.694, 43.9, 1.7]
# μ, b, θ, cₓ, cₚ, ν, w₀, γ, f₀, c_f, β = p0


ec = ContinuousDynamicalSystem(excitable_carbon, u0, p0)

# streamplot:
xgrid = 0:5:190
ygrid = 2000:20:3500
ux = zeros(length(xgrid), length(ygrid))
uy = copy(ux)
for (i, x) in enumerate(xgrid)
    for (j, y) in enumerate(ygrid)
        ux[i, j], uy[i, j] = ec.f(SVector(x, y), ec.p, 0)
    end
end

fig = figure()
ax = fig.add_subplot(1,3,1)
ax.streamplot(Vector(xgrid), Vector(ygrid) ./1000 , ux', uy' ./ 1000;
    linewidth = 1.0, density = 0.25, color = "C5", arrowsize = 1.0
)
# random trajectory
c_stable = 83.58 # approximate value
for (i, u) in enumerate(([c_stable, 2075.0], [c_stable, 2102.0], ))
    tr = trajectory(ec, 100.0, u)
    c, w = columns(tr)
    w ./= 1000
    ax.plot(c, w; color = "C$i")
end
ax.set_xlabel("\$c\$"; labelpad = -15)
ax.set_ylabel("\$w \\times 10^{-3}\$")
ax.set_xticks(0:60:180)

fig.tight_layout(pad = 0.3)
add_identifiers!(fig)
wsave(plotsdir("carboncycle"), fig)



# %% Fractal dimension of chaotic attractors
fig = figure()
ax_ts = fig.add_subplot(1,3,1)
ax_3d = fig.add_subplot(1,3,2; projection = "3d")
ax_fd = fig.add_subplot(1,3,3)

vostok = readdlm(projectdir("exercise_data", "17.csv"))
time = vostok[:, 1]
δT = vostok[:, 2]

# ax_ts.plot(time, δT; lw = 1.0, alpha = 0.5)

# Interpolate Vostok data
using Dierckx
# interpolation object:
spl = Spline1D(time, δT; k=3, bc="extrapolate")
# equi-spaced time vector:
t = range(minimum(time), maximum(time); length = 2000)

T = spl(t)

ax_ts.plot(t ./ 1000, T; lw = 1.0)
ax_ts.set_xlabel("time (towards past, kyr)")
ax_ts.set_ylabel("\$\\delta T\$ (K)")

# mi = selfmutualinfo(T, 1:1000)
# ax_fd.plot(1:1000, mi)
τ = 75 # inspect the above plot to see a good delay time

# plot attractor
R = embed(T, 3, τ)
ax_3d.plot(columns(R)...; lw = 1.0)
ax_3d.set_xticklabels([])
ax_3d.set_yticklabels([])
ax_3d.set_zticklabels([])

# estimate fractal dim
ds = 3:12
colors = matplotlib.cm.inferno(ds ./ (ds[end] + 2))
es = estimate_boxsizes(R; z = 0)

for (i, d) in enumerate(ds)
    A = embed(T, d, τ)
    Cs = boxed_correlationsum(A, es)
    j = findfirst(z -> z > 0, Cs)
    x, y = log.(es)[j:end], log.(Cs)[j:end]
    # x, y = log.(es), log.(Cs)
    # li, Δ = linear_region(x, y)
    c = colors[i, :]
    ax_fd.plot(x, y; color = c, lw = 1.5)
    # li = findall(z -> -4 ≤ z ≤ 0, x)
    _, Δ = linreg(x, y)
    @show (d, Δ)
    # ax_fd.plot(x[[1, end]], y[[1, end]],
    # zorder = 2, ms = 10, marker = "o", ls = "None", color = c, alpha = 0.75,
    # label = "d=$(d), Δ=$(round(Δ; digits=2))", )
end

# Show two regions
# ax_fd.axvspan(-3, -0.25; color = "C1", alpha = 0.1)
el1 = matplotlib.patches.Ellipse(
    (-1, -10), 3, 5; angle = -15, linewidth = 2, color = "C1",
    fill = false
)
ax_fd.add_patch(el1)
el2 = matplotlib.patches.Ellipse(
    (0.5, -5), 4, 6; angle = -20, linewidth = 2, color = "C2",
    fill = false
)
ax_fd.add_patch(el2)

ax_fd.set_xlabel("\$\\log(\\varepsilon)\$")
ax_fd.set_ylabel("\$\\log(C)\$")
add_identifiers!(fig)
fig.tight_layout(;pad = 0.3)
fig.tight_layout(;pad = 0.3)
# wsave(plotsdir("vostok"), fig)
