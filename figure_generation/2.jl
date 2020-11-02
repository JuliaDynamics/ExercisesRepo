using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# %% 1d climate
using DynamicalSystems, PyPlot, Roots

αtan(T) = 0.5 - (0.4/π)*atan((T-263)/2)
dTdt(T, ε = 0.65, α=αtan, s= 1.0) = s*(1 - α(T)) - 1.6e-10 * ε * T^4
dTdt(T; ε = 0.65, α=αtan, s = 1.0) = dTdt(T, ε, α, s)

# Ts = 200:400.0
# plot(Ts, dTdt.(Ts))
# plot(Ts, dTdt.(Ts, 0.2))
# plot(Ts, dTdt.(Ts, 0.9))
# axhline(0)

close("all")
fig = figure(figsize = (figx/1.5, figx/3))
Ts = 200:0.5:320.0
arrows = 210:10:300 |> collect
deleteat!(arrows, findfirst(isequal(260), arrows))
roots = Roots.find_zeros(dTdt, Ts[1], Ts[end])
plot(Ts, dTdt.(Ts), color = "C1", label = "\$dT/dt\$")
axhline(0; lw = 2.0, zorder = -99)
xlim(Ts[1], Ts[end])
ylim(-0.2, 0.2)
yticks([-0.1, 0, 0.1])
xlabel("\$T\$")
tight_layout()
for (i, r) in enumerate(roots)
    plot(r, 0, marker = "o", markeredgecolor = "k", markerfacecolor = iseven(i) ? "w" : "k",
    markersize = 15, mew = 2)
end
for r in arrows
    f = dTdt(r)
    x, dx = f > 0 ? (r - 5, 5) : (r+5, -5)
    ff = abs(1.2f)^2 + 0.01
    arrow(x, 0, dx, 0; color = "C2", width = ff, length_includes_head = false,
    head_width = 1.5ff, head_length = 100ff
    )
end
legend()
tight_layout()
subplots_adjust(bottom = 0.2, left = 0.1)
# fsave(joinpath(figdir, "1dstatespace"), fig)

# %% Bifurcation versus ε
figure()
es = 0.2:0.01:1.0
for e in es
    f = (T) -> dTdt(T, e)
    roots = Roots.find_zeros(f, Ts[1], Ts[end])
    for (i, r) in enumerate(roots)
        plot(e, r, marker = "o", markeredgecolor = "k", markerfacecolor = iseven(i) ? "w" : "k",
        markersize = 10, mew = 2)
    end
end
plot([],[], ls = "None", marker = "o", markeredgecolor = "k", markerfacecolor = "k", label = "stable")
plot([],[], ls = "None", marker = "o", markeredgecolor = "k", markerfacecolor = "w", label = "unstable")
legend()
xlabel("\$\\epsilon\$")
ylabel("\$T^*\$")
tight_layout()

# %% FitzHugh
fig, axs = subplots(1,3; sharey = true)
Is = [0.2, 0.2, 0.5]
ds = Systems.fitzhugh_nagumo(; I = Is[1])
xgrid = -2.5:0.1:2.5
ygrid = -1.0:0.1:1.8

function add_nullclines!(ax, I)
    v = -2.5:0.1:2.5
    w1 = @. (v + 0.7)/0.8
    w2 = @. v - v^3/3 + I
    ax.plot(v, w1, color = "C3")
    ax.plot(v, w2, color = "C2")
    ax.set_ylim(ygrid[1], ygrid[end])
end

function add_streamlines!(ax, I)
    ux = zeros(length(xgrid), length(ygrid))
    uy = copy(ux)
    set_parameter!(ds, 1, I)
    for (i, x) in enumerate(xgrid)
        for (j, y) in enumerate(ygrid)
            ux[i, j], uy[i, j] = ds.f(SVector(x, y), ds.p, 0)
        end
    end
    ax.streamplot(Vector(xgrid), Vector(ygrid), ux', uy'; linewidth = 0.5, color = "C0")
    ax.set_ylim(ygrid[1], ygrid[end])
end

for i in 1:3;
    axs[i].set_title("\$I = $(Is[i])\$");
    i == 1 && add_streamlines!(axs[i], Is[i])
    add_nullclines!(axs[i], Is[i])
end

# Add tex to axs [1]
sss = 30
props = Dict(:boxstyle=>"round", :alpha=>0.75, :facecolor => "white")

axs[1].text(0.9, 1.4, "\$\\dot{w}>0\$", color = "C3", size = sss, bbox=props)
axs[1].text(0.1, 1.4, "\$\\dot{w}<0\$", color = "C3", size = sss, ha = :right, bbox=props)
axs[1].text(-1.5, 0.7, "\$\\dot{v}<0\$", color = "C2", size = sss, bbox=props)
axs[1].text(-0.3, -0.6, "\$\\dot{v}>0\$", color = "C2", size = sss, bbox=props)

for i in 2:3;
    u1 = [-1.8, -0.5]
    u2 = [-0.45, -0.6]
    u3 = [-0.15, -0.2]
    for (j, u) in enumerate((u1, u2))
        set_parameter!(ds, 1, Is[i])
        tr = trajectory(ds, 100.0, u)
        axs[i].plot(columns(tr)...; ls = j == 1 ? "-" : "--", color = "C$(j-1)")
        axs[i].scatter(tr[1]...; color = "C$(j-1)", s = 80)
    end
end

axs[1].set_ylabel("\$w\$")
axs[2].set_xlabel("\$v\$")
# axs[3].xaxis.set_label_coords(1.0, 0.0)

fig.tight_layout()
fig.subplots_adjust(bottom = 0.18, left = 0.08, wspace = 0.1)
add_identifiers!(fig)
fsave(joinpath(figdir, "fitzhugh"), fig)

# %% Hose noover
using DynamicalSystems, PyPlot
# Generate PSOS first, to be able to pick points
ds = Systems.nosehoover()
psos = poincaresos(ds, (1, 0.0), 10000.0)
scatter(psos[:, 2], psos[:, 3])

ch = [0, 0.1, 0]
qp = [0, 1.0, 0]
pd =  [0, 1.549934227822898, 0]

psos = poincaresos(ds, (1, 0.0), 10000.0; qp)
scatter(psos[:, 2], psos[:, 3])

# Okay now plot the actual 3D stuff
# %%
using3D()
figure()

# using Makie
# sc = Scene()
for (i, u0) in enumerate((qp, pd))

    tr = trajectory(ds, 1000, u0; dt = 0.05)
    plot3D(columns(tr)...; linewidth = i^2*1.0,)
    # Makie.lines!(sc, columns(tr)...; linewidth = i^3*1.0, color = to_color(COLORS[i]))
end
# display(sc)

# %% Power spectra for three orbits of Henon-Heiles
using DynamicalSystems, PyPlot, FFTW, Statistics

u0s = (
    [0.0, -0.25, 0.42, 0.0],
    [0.0, 0.30266571044921875, 0.4205654433900762, 0.0],
    [0.0, 0.1, 0.5, 0.0],
)

labels = (
    "chaotic",
    "periodic",
    "quasiperiodic",
)

hh = Systems.henonheiles()
fig, axs = subplots(3, 1; figsize = (figx/2, 2figy), sharex = true)

δt = 0.05

for (i, u) in enumerate(u0s)
   # r = trajectory(hh, 1000.0, u; dt = 0.1)[:, 1]
    r = trajectory(hh, 30000.0, u; dt = δt)[:, 1]
    P = abs2.(rfft(r .- mean(r)))
    P[1] = P[2]
    ν = rfftfreq(length(r))/δt
    # axs[i].plot(ν, P ./ maximum(P),
    axs[i].semilogy(ν, P ./ maximum(P),
    label = labels[i], linewidth = i == 1 ? 1.0 : 2.0, color = "C$(i-1)")
    axs[i].text(0.99, 0.8, labels[i]; ha = "right", transform = axs[i].transAxes,
    color = "C$(i-1)")
    # @show std(r)
    # r .+= 0.5randn(length(r))
    # P = abs2.(rfft(r .- mean(r)))
    # PyPlot.plot(10rfftfreq(length(r)), P ./ maximum(P),
    # lw = 1.0, alpha = 0.5, color = "C$(i-1)")
    # PyPlot.plot(r)
    # axs[i].set_yticks([])
    axs[i].set_ylim(10.0^(-5), 1.0)
    axs[i].set_yticks(10.0 .^ (-1:-2:-5))
    axs[i].set_xlim(0, 0.4)
end
axs[2].set_ylabel("\$P / \\mathrm{max}(P)\$")
axs[3].set_xlabel("frequency \$\\nu\$")
fig.tight_layout()
fig.subplots_adjust(bottom = 0.1, left = 0.2, top = 0.98, right = 0.93, hspace = 0.1)
fsave(plotsdir("spectra"), fig)

# %% fitzhugh new
fig, axs = subplots(2,3; figsize = (figx, 2figy))
fig.tight_layout()

a = 3.
b = 0.2
ε = 0.01
I = 0.0

# ds = Systems.fitzhugh_nagumo_new(; a, b, ε, I )
ds = Systems.fitzhugh_nagumo([0,0.0]; a, b, ε, I )

xgrid = -0.3:0.02:1.1
ygrid = -0.1:0.01:0.4

function add_nullclines!(ax, a, b, ε, I)
    u = -0.3:0.01:1.1
    w1 = @. u
    w2 = @. a*u*(u-b)*(1-u) + I
    ax.plot(u, w1, "--", color = "C3")
    ax.plot(u, w2, ":", color = "C4")
    ax.set_ylim(ygrid[1], ygrid[end])    # why is this necessary?
end

function add_streamlines!(ax, a, b, ε, I)
    ux = zeros(length(xgrid), length(ygrid))
    uy = copy(ux)
    set_parameter!(ds, [a, b, ε, I])
    for (i, x) in enumerate(xgrid)
        for (j, y) in enumerate(ygrid)
            ux[i, j], uy[i, j] = ds.f(SVector(x, y), ds.p, 0)
        end
    end
    ax.streamplot(Vector(xgrid), Vector(ygrid), ux', uy'; linewidth = 0.5, density = 0.5, color = "C7")
    ax.set_ylim(ygrid[1], ygrid[end])
end

# Fig A (top left)
# ================
ax = axs[1]
add_nullclines!(ax, a, b, ε, I)
add_streamlines!(ax, a, b, ε, I)

# Add tex to axs [1] explaining nullclines
sss = 20
props = Dict(:boxstyle=>"round", :alpha=>1., :facecolor => "white")
ax.text(0.5, 0.35, "\$\\dot{w}>0\$", color = "C3", size = sss, bbox=props)
ax.text(-0.1, 0.35, "\$\\dot{w}<0\$", color = "C3", size = sss, bbox=props)
ax.text(0, 0.1, "\$\\dot{u}<0\$", color = "C4", size = sss, bbox=props)
ax.text(0.45, -0.05, "\$\\dot{u}>0\$", color = "C4", size = sss, bbox=props)
ax.set_ylabel("\$w\$")
ax.set_xlabel("\$u\$")


# Fig B & C: two trajectries & zoomin
# ==================
for ax in (axs[3], axs[5])
    add_nullclines!(ax, a, b, ε, I)

    u1 = [0.19, 0.]
    tr1 = trajectory(ds, 200.0, u1)
    ax.plot(columns(tr1)...; color = "C2", lw = 1.5)
    ax.scatter(tr1[1]...; color = "C2", s = 20)

    u2 = [0.21, 0.]
    tr2 = trajectory(ds, 300.0, u2)
    ax.plot(columns(tr2)...; color = "C1", lw = 1.5)
    ax.scatter(tr2[1]...; color = "C1", s = 20)

    # ic = (tr1[1] + tr2[1]) ./ 2
    # ax.scatter(ic...; color = "k", s = 10)

    ax.plot(0, 0; marker = "o", mec = "C0", mew = 1,
        markersize = 6, mfc = "w", zorder = 99)
    ax.axhline(0; color = "k", lw = 0.5, zorder = 1)
end
axs[3].set_ylabel("\$w\$")
axs[3].set_xlabel("\$u\$")

zbox = ((-0.05, -0.04), (0.25, 0.04))
axis_zoomin!(axs[5],axs[3],  zbox, zbox, "C0"; lw = 1.0)
axs[5].set_xlim(zbox[1][1], zbox[2][1])
axs[5].set_ylim(zbox[1][2], zbox[2][2])
axs[5].set_xticks([])
axs[5].set_yticks([])

# Fig D: timeseries, pulses, refractory period
using OrdinaryDiffEq
pulses_start = [20, 80, 170] # callback times
pulse_width = 4
pulses_end = pulses_start .+ pulse_width
pulses = sort!(vcat(pulses_start, pulses_end))
Ipulse = 0.2
condition(u,t,integ) = t ∈ pulses # pulse times
function affect!(integ)
    i = integ.t ∈ pulses_start ? Ipulse : 0.0
    integ.p[4] = i
end
cb = DiscreteCallback(condition, affect!)

Tf = 250.0
dt = 0.1
prob = ODEProblem(ds, (0.0, Tf))
sol = solve(prob, Tsit5(); callback=cb, tstops = pulses, dtmax = 0.01, maxiter = typemax(Int))

axs[2].plot(sol.t, sol[1, :])
pt = [any(x -> x ≤ t ≤ x + pulse_width, pulses_start) ? Ipulse : 0.0 for t in sol.t]
axs[2].plot(sol.t, pt; color = "C2")
axs[2].set_xlabel("\$t\$")
axs[2].set_ylabel("\$u\$")

# Fig E: multi stable
a = 8.
b = 0.2
ε = 0.01
I = 0.0

xgrid = -0.3:0.02:1.1
ygrid = -0.2:0.01:1

add_nullclines!(axs[4], a, b, ε, I)
set_parameter!(ds, [a,b,ε,I])
ax = axs[4]

u1 = [0.3, 0.]
tr1 = trajectory(ds, 300.0, u1)
ax.plot(columns(tr1)...; color = "C2")
ax.scatter(tr1[1]...; color = "C2", s = 20)

u2 = [0.55, 0.8]
tr2 = trajectory(ds, 300.0, u2)
ax.plot(columns(tr2)...; color = "C1")
ax.scatter(tr2[1]...; color = "C1", s = 20)

ax.set_ylabel("\$w\$")
ax.set_xlabel("\$u\$")

# Fig F: limit cycle
# -------------------------
ax = axs[6]
a = 3.
b = -0.05
ε = 0.01
I = 0.
set_parameter!(ds, [a,b,ε,I])
add_nullclines!(ax, a, b, ε, I)

tr = trajectory(ds, 400.0, [0.0, 0.1])
ax.plot(columns(tr)...; color = "C2")
ax.scatter(tr[1]...; color = "C2", s = 20)
ax.set_xlim(-0.5,1.1)
ax.set_ylim(-0.1,0.6)
ax.set_ylabel("\$w\$")
ax.set_xlabel("\$u\$")
