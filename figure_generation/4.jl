using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))

# %% 1d climate
using DynamicalSystems, PyPlot, Roots, ForwardDiff

# %% Bifurcation versus ε
fig = figure(figsize = (figx/2, figx/3))
ax = gca()

α(T) = 0.5 - 0.2*tanh((T-263)/4)
dTdt(T, ε=0.65) = 1 - α(T) - 10*ε*(0.002T)^4
es = 0.3:0.005:0.9
for e in es
    f = (T) -> dTdt(T, e)
    roots = Roots.find_zeros(f, 50, 500)
    for (i, r) in enumerate(roots)
        c = ForwardDiff.derivative(f, r) < 0 ? "k" : "w"
        plot(e, r, marker = "o", mec = "k", mew = 1,
        markersize = 6, mfc = c )
    end
end
plot([],[], ls = "None", marker = "o", mec = "k", mfc = "k", label = "stable")
plot([],[], ls = "None", marker = "o", mec = "k", mfc = "w", label = "unstable")
legend(markerscale = 3)
xlabel("\$\\epsilon\$")
ylabel("\$T^*\$")
tight_layout()

e1, T1 = 0.48755, 256.38
e2, T2 = 0.78661, 268.075
for (e, T) in ((e1, T1), (e2, T2))
    # plot(e, T; marker="o", ms = 10, color = "k", mew = 1, mec = "C2",
    #            fillstyle = "top", markerfacecoloralt= "w")
    plot(e, T; marker="o", ms = 10, color = "C2")
end

fig.subplots_adjust(left = 0.18, bottom = 0.2, top = 0.97)
# fsave(joinpath(figdir, "bif_example"), fig)

# %% Period doubling

# %% Roessler trajectory -> PSOS -> Lorenz map
using InteractiveDynamics
using DynamicalSystems, GLMakie
using AbstractPlotting.MakieLayout


N = 10000.0

ds = Systems.roessler()
tr = trajectory(ds, N; Ttr = 100.0)

figure, layout = layoutscene(resolution = (2100, 700), )
display(figure)

# Plot trajectory

trplot = layout[1,1] = LScene(figure, scenekw = (camera = cam3d!, raw = false))
lines!(trplot, columns(tr)...; color = COLORSCHEME[1], linewidth = 2.0)

# Adjust ticks and sizes of 3D plot
xticks!(trplot.figure; xtickrange = [-10, 0, 10], xticklabels = string.([-10, 0, 10]))
yticks!(trplot.figure; ytickrange = [-10, 0, 10], yticklabels = string.([-10, 0, 10]))
zticks!(trplot.figure; ztickrange = [5, 15, 25], zticklabels = string.([5, 15, 25]))
trplot.figure[Axis][:ticks][:textsize] = (15, 15, 15)
trplot.figure[Axis][:names][:textsize] = (20,20,20)

# Plot plane and section
o = Point3f0(-10, 0, 0)
w = Point3f0(25, 0, 25)
p = FRect3D(o, w)
a = RGBAf0(0,0,0,0)
c = to_color(COLORSCHEME[2])
img = AbstractPlotting.ImagePattern([c a; a c]); # This throws an error if it shows, you can ignore that
mesh!(trplot, p; color = img)

psos = poincaresos(ds, (2, 0.0), N; Ttr = 100.0)
scatter!(trplot, columns(psos)...; color = COLORSCHEME[3], markersize = 500)

# Plot section separately
psplot = layout[1, 2] = LAxis(figure)
scatter!(psplot, psos[:, 1], psos[:, 3]; color = COLORSCHEME[3])
psplot.xlabel = "xₙ"
psplot.ylabel = "zₙ"

LS = 50
TS = 40

psplot.xticklabelsize = TS
psplot.yticklabelsize = TS
psplot.xlabelsize = LS
psplot.ylabelsize = LS

# Plot lorenz map
loplot = layout[1, 3] = LAxis(figure)
scatter!(loplot, psos[1:end-1, 3], psos[2:end, 3]; color = COLORSCHEME[3])
loplot.xlabel = "zₙ"
loplot.ylabel = "zₙ₊₁"
loplot.xticklabelsize = TS
loplot.yticklabelsize = TS
loplot.xlabelsize = LS
loplot.ylabelsize = LS
display(figure)

# %% Bifurcation kit code
# Better check https://rveltz.github.io/BifurcationKit.jl/dev/iterator/#
using BifurcationKit, SparseArrays, LinearAlgebra
const BK = BifurcationKit
using ForwardDiff

F(T, ε=0.65) = @. 0.5 + 0.2*tanh((T-263)/4) - 10*ε*(0.002T)^4

function deriv(x, p)
	f = (x) -> F(x, p)
	return ForwardDiff.derivative.(Ref(f), x)
end
Jac_m = (x, p) -> diagm(0 => deriv(x, p))


# parameters for the continuation
opts = ContinuationPar(
	dsmax = 0.1, dsmin = 1e-3,
	ds = 0.01, maxSteps = 2000, pMin = 0.0, pMax = 0.9, saveSolEveryStep = 0,
	newtonOptions = BK.NewtonPar(tol = 1e-8, verbose = false),
	detectBifurcation = true
)

# we define an iterator to hold the continuation routine
p0 = 0.3
u0 = [350.0]
iter = BK.PALCIterable(F, Jac_m, u0, p0, (BK.@lens _), opts; verbosity = 2)

resp = Float64[]
resx = Float64[]
lambdas = Float64[]

# this is the PALC algorithm
for state in iter
	# we save the current solution on the branch
	push!(resx, getx(state)[1])
	push!(resp, getp(state))
	BK.computeEigenvalues!(iter, state)
	eig = state.eigvals
	push!(lambdas, real(eig[1]))
end

fig = figure(figsize = (figx/2, figx/3))
ax = gca()
i1 = findfirst(l -> l ≥ 0, lambdas)
i2 = findlast(l -> l ≥ 0, lambdas)
ax.plot(resp[1:i1], resx[1:i1], color = "C3", label = "stable")
ax.plot(resp[i1:i2], resx[i1:i2], color = "C3", ls = "--", label = "unstable")
ax.plot(resp[i2+1:end], resx[i2+1:end], color = "C3")
ax.plot([resp[i1], resp[i2]], [resx[i1], resx[i2]], ls = "None",
		marker="o", ms = 10, color = "k", mew = 1, mec = "k",
        fillstyle = "top", markerfacecoloralt= "w"
)

# PyPlot.plot(resp, resx)
ax.legend(markerscale = 3)
ax.set_xlabel("\$\\epsilon\$")
ax.set_ylabel("\$T^*\$")
fig.tight_layout()
fig.subplots_adjust(left = 0.18, bottom = 0.2, top = 0.97)
# fsave(joinpath(figdir, "bif_example"), fig)


# %% Intermittency simple Manneville map
using DynamicalSystems, PyPlot
lo = Systems.manneville_simple(0.454)

rs = [-0.1, 0.0001, 0.05]

fig = figure(figsize = (figx, 1.5figy))
ax1 = fig.add_subplot(3, 2, 1)
ax2 = fig.add_subplot(3, 2, 3; sharex = ax1)
ax3 = fig.add_subplot(3, 2, 5; sharex = ax1)
axs = (ax1, ax2, ax3)
axc = fig.add_subplot(1, 2, 2)
add_identifiers!(fig)

T = 180
xs = 0:0.0001:1

for i in 1:length(rs)
	set_parameter!(lo, 1, rs[i])
	x = trajectory(lo, T; Ttr = 100)
	chaotici = findall(a -> a>0.1, x)
	regulari = setdiff(1:length(x), chaotici)
	axs[i].plot(0:length(x)-1, x; lw = 0.5, color = "C0")
	axs[i].plot(regulari .- 1, x[regulari]; ls = "None", marker = "o", color = "C0")
	axs[i].plot(chaotici .- 1, x[chaotici]; ls = "None", marker = "o", color = "C2")
	axs[i].set_ylim(-0.1, 1.2)
	axs[i].set_yticks([0, 1])
	# axs[i].text(0.99, 0.90, "\$r=$(rs[i])\$"; bbox = bbox,
	# transform = axs[i].transAxes, va = :top, ha=:right, fontsize = 22)
end
# axs[1].set_xlim(0, T)
for ax in (ax1, ax2); ax.tick_params(labelbottom=false); end
axs[1].set_xlim(0, T)
axs[1].set_xticks(0:60:180)
axs[3].set_xlabel("\$n\$"; labelpad = -25)
axs[2].set_ylabel("\$x_n\$"; labelpad = -10)
axs[2].axvspan(21, 40; color = "C3", alpha = 0.5)

# laminar arrow
xc = (135 + 169)/2
xspan = (169 - 135)
nice_arrow!(ax2, xc, 1.0, xspan, 0; tex = "\$\\ell\$", yo = 0.0, xo = -xspan/2 - 5)

xc = (65 + 100)/2
xspan = abs(65 - 100)
nice_arrow!(ax2, xc, 1.0, xspan, 0; tex = "\$c\$", yo = 0.0, xo = xspan/2, color = "C2")

# Do the cobweb plot
function cobweb(t) # transform timeseries t into cobweb (point2D)
    cx = Float64[]; cy = Float64[]
    for i ∈ 1:length(t)-1
        push!(cx, t[i]); push!(cy, t[i])
		push!(cx, t[i]); push!(cy, t[i+1])
    end
    return cx, cy
end


axc.set_xlim(0.0, 1.0)
axc.set_ylim(0.0, 1.0)
axc.set_xticks([0,1])
axc.set_yticks([0,1])
axc.set_ylabel("\$x_{n+1}\$"; labelpad = -15)
axc.set_xlabel("\$x_{n}\$"; labelpad = -25)
axi = axc.inset_axes([0.05, 0.55, 0.3, 0.4])

for ax in (axc, axi)
	_r = rs[2]
	set_parameter!(lo, 1, _r)
	f = lo.f.(xs, Ref([_r]), 0)
	ax.plot(xs, f, color = "C1")
	ax.plot([0, 1], [0,1]; lw = 1.0, color = "k")
	x = trajectory(lo, 60; Ttr = 100)
	cx, cy = cobweb(x)
	ax.plot(cx, cy; lw = 1.0, alpha = 0.5, color = COLORS[1])
end

z = 0.01
w = 0.05
zbox = ((z, z), (z+w, z+w))
axis_zoomin!(axi, axc, zbox, zbox, "C3"; dir = :top)
axi.set_ylim(z, z+w)
axi.set_xlim(z, z+w)
axi.set_xticks([])
axi.set_yticks([])
axc.grid(false)

fig.tight_layout(pad = 0.3)
# wsave(plotsdir("intermittency_manneville"), fig)
