using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))

# %% 1d climate
using DynamicalSystems, PyPlot, Roots, ForwardDiff

# %% Bifurcation versus ε
fig = figure(figsize = (figx/2, figx/3))
ax = gca()
dTdt(T, ε=0.65) = 0.5 + (0.4/π)*atan((T-263)/2) - 10*ε*(0.002T)^4
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
using DynamicalSystems, PyPlot

fig = figure(figsize = (figx, 1.2figx))
axlo = fig.add_subplot(3, 1, 1)
axro = fig.add_subplot(3, 1, 3)

zoomaxs = []
for i in 1:3
	push!(zoomaxs, fig.add_subplot(3, 3, 3+i))
end

axlo.clear()
ds = Systems.logistic(0.5)
i = 1
pvalues = range(2.5, 4; length = 2001)
n = 2000
Ttr = 2000
p_index = 1
output = orbitdiagram(ds, i, p_index, pvalues; n = n, Ttr = Ttr)

L = length(pvalues)
x = Vector{Float64}(undef, n*L)
y = copy(x)
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= pvalues[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end
axlo.plot(x, y, ls = "None", ms = 0.5, color = "black", marker = "o", alpha = 0.01)
axlo.set_xlim(pvalues[1], pvalues[end]); axlo.set_ylim(0,1)
axlo.set_ylabel("\$x_\\infty\$", labelpad = 10);
axlo.text(2.6, 0.1, "logistic map", size = 36)
r∞ = 3.56995
axlo.axvline(r∞; color = "C3", alpha = 0.7, ls = "dashed")

# Zoomin boxes
zbox1 = ((3.629, 0.295), (3.640, 0.32))
zbox2 = ((3.6338, 0.304), (3.6342, 0.3065))
zbox3 = ((3.63407, 0.30485), (3.63408, 0.3058))
allboxes = (zbox1, zbox2, zbox3)
for (k, zbox) in enumerate(allboxes)
	pval = range(zbox[1][1], zbox[2][1]; length = 1001)
	nsm = n÷3
	output = orbitdiagram(ds, i, p_index, pval;
		n=nsm, Ttr, ulims = (zbox[1][2], zbox[2][2])
	)
	L = length(pval)
	x = Vector{Float64}(undef, nsm*L)
	y = copy(x)
	for j in 1:L
	    x[(1 + (j-1)*nsm):j*nsm] .= pval[j]
	    y[(1 + (j-1)*nsm):j*nsm] .= output[j]
	end
	zoomaxs[k].clear()
	zoomaxs[k].plot(x, y, ls = "None", ms = 0.5, color = "black", marker = "o", alpha = 0.01)
	zoomaxs[k].set_xlim(zbox[1][1], zbox[2][1])
	zoomaxs[k].set_ylim(zbox[1][2], zbox[2][2])
	zoomaxs[k].axis("off")
end

zbox0 = ((pvalues[1], 0.0), (pvalues[end], 1.0))
allzoomaxs = [axlo, zoomaxs...]
for i in 1:3
	origin = allzoomaxs[i]
	zoomin = allzoomaxs[i+1]
    zbox = allboxes[i]
	c = i != 3 ? "C$(i)" : "C$(i+2)"
    axis_zoomin!(zoomin, origin, zbox, zbox, c, 1.0,
	connect_lines = i > 1, lw = 4.0)
end

# arrow indicator
axlo.arrow(
	(zbox1[1][1] + zbox1[2][1])/2, 0.05, 0, zbox1[2][2] - 0.1;
	color = "C1", width = 0.005, length_includes_head = true, head_length = 0.1
)

# Roessler
ro = Systems.roessler()
pvalues = range(1.5, stop = 5.18, length = 2001)
i = 1
plane = (2, 0.0)
tf = 4000.0
p_index = 3
output = produce_orbitdiagram(ro, plane, i, p_index, pvalues;
                              tfinal = tf, Ttr = 2000.0)

for (j, p) in enumerate(pvalues)
    axro.plot(fill(p, length(output[j])), output[j], lw = 0,
    marker = "o", ms = 0.2, color = "black", alpha = 0.1)
end
axro.set_xlim(pvalues[1], pvalues[end]);
axro.set_ylim(minimum(minimum(o) for o in output), maximum(maximum(o) for o in output))
axro.set_ylabel("\$x_\\infty\$ (section)", labelpad = 25)
axro.set_xlabel("parameter (\$r\$ for logistic, \$c\$ for Rossler)")
axro.set_yticks(3:3:9)
axro.text(1.7, 7, "Rössler system", size = 36)
fig.tight_layout()
fig.subplots_adjust(bottom = 0.1, top = 0.97, left = 0.09, hspace = 0.3, wspace = 0.1, right = 0.98)
# fsave(plotsdir("orbit_diagrams"), fig)

# %% # Testing how orbit diagram of quasiperiodic looks like
using DynamicalSystems, PyPlot

hh = Systems.henonheiles([0.0, 0.1, 0.5, 0.0])
psos = poincaresos(hh, (1, 0.0), 5000.0)

figure()
for i in 0:0.001:1
    plot(fill(i, length(psos)), psos[:, 2], lw = 0,
    marker = "o", ms = 0.2, color = "black", alpha = 0.1)
end

sm = Systems.standardmap([π/2, 0.1])
n = 1000
Ttr = 2000
pvalues = range(0, 1.2; length = 2001)
output = orbitdiagram(sm, 2, 1, pvalues; n = n, Ttr = Ttr)
figure()
ax = gca()
L = length(pvalues)
x = Vector{Float64}(undef, n*L)
y = copy(x)
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= pvalues[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end
ax.plot(x, y, ls = "None", ms = 0.5, color = "black", marker = "o", alpha = 0.02)

# %% Roessler trajectory -> PSOS -> Lorenz map
using InteractiveChaos
using DynamicalSystems, Makie
using AbstractPlotting.MakieLayout


N = 10000.0

ds = Systems.roessler()
tr = trajectory(ds, N; Ttr = 100.0)

scene, layout = layoutscene(resolution = (2100, 700), )
display(scene)

# Plot trajectory

trplot = layout[1,1] = LScene(scene, scenekw = (camera = cam3d!, raw = false))
lines!(trplot, columns(tr)...; color = COLORSCHEME[1], linewidth = 2.0)

# Adjust ticks and sizes of 3D plot
xticks!(trplot.scene; xtickrange = [-10, 0, 10], xticklabels = string.([-10, 0, 10]))
yticks!(trplot.scene; ytickrange = [-10, 0, 10], yticklabels = string.([-10, 0, 10]))
zticks!(trplot.scene; ztickrange = [5, 15, 25], zticklabels = string.([5, 15, 25]))
trplot.scene[Axis][:ticks][:textsize] = (15, 15, 15)
trplot.scene[Axis][:names][:textsize] = (20,20,20)

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
psplot = layout[1, 2] = LAxis(scene)
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
loplot = layout[1, 3] = LAxis(scene)
scatter!(loplot, psos[1:end-1, 3], psos[2:end, 3]; color = COLORSCHEME[3])
loplot.xlabel = "zₙ"
loplot.ylabel = "zₙ₊₁"
loplot.xticklabelsize = TS
loplot.yticklabelsize = TS
loplot.xlabelsize = LS
loplot.ylabelsize = LS
display(scene)

# Makie.save(plotsdir("roessler_map.png"), scene)

# %% Bifurcation kit code
# Better check https://rveltz.github.io/BifurcationKit.jl/dev/iterator/#
using BifurcationKit, SparseArrays, LinearAlgebra
const BK = BifurcationKit
using ForwardDiff

F(T, ε=0.65) = @. 0.5 + (0.4/π)*atan((T-263)/2) - 10*ε*(0.002T)^4

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

# %% Intermittency in the Roessler system
using DynamicalSystems, PyPlot
ro = Systems.roessler(ones(3))
cc = 5.187

cs = range(5.187; step = -0.002, length = 3)

fig, axs = subplots(length(cs), 1; sharex = true)
T = 3000
dt = 0.05
for i in 1:length(cs)
	set_parameter!(ro, 3, cs[i])
	tr = trajectory(ro, T; Ttr = 200.0, dt)
	axs[i].plot(0:dt:T, tr[:, 1], lw = 0.5)
	axs[i].set_yticks([])
	axs[i].text(0.99, 0.90, "\$c=$(cs[i])\$"; bbox = bbox,
	transform = axs[i].transAxes, va = :top, ha=:right, fontsize = 22)
end
axs[1].set_xlim(0, T)
fig.tight_layout(pad = 0.2)
fig.subplots_adjust(hspace = 0.1, bottom = 0.1)

# %% Intermittency in the PSOS of Roessler system
using DynamicalSystems, PyPlot

cc = 5.187
plane = (2, 0.0)

ro = Systems.roessler(ones(3))
cs = range(5.187; step = -0.002, length = 3)

fig, axs = subplots(length(cs), 1; sharex = true)
T = 3000
dt = 0.05
for i in 1:length(cs)
	set_parameter!(ro, 3, cs[i])
	psos = poincaresos(ro, plane, T; Ttr = 200.0, idxs = [1])
	axs[i].plot(psos[:, 1]; lw = 0.5)
	axs[i].set_yticks([])
	axs[i].text(0.99, 0.90, "\$c=$(cs[i])\$"; bbox = bbox,
	transform = axs[i].transAxes, va = :top, ha=:right, fontsize = 22)
end
# axs[1].set_xlim(0, T)
fig.tight_layout(pad = 0.2)
fig.subplots_adjust(hspace = 0.1, bottom = 0.1)

# %% Intermittency in logistic map
using DynamicalSystems, PyPlot
lo = Systems.logistic(0.4)
r3 = 1 + sqrt(8)
rs = [r3, r3 - 0.0001, r3 - 0.001]
# rs = range(r3; step = -0.0001, length = 3)

chaoticphases = (
[],
[27:42, 141:152, 170:200],
[0:6, 33:50, 55:76, 83:89, 113:145, 161:200]
)

fig, axs = subplots(length(rs), 1; sharex = true)
T = 200
dt = 0.05
for i in 1:length(rs)
	set_parameter!(lo, 1, rs[i])
	x = trajectory(lo, T; Ttr = 100)
	axs[i].plot(x; lw = 0.5, marker = "o")
	for idx in chaoticphases[i]
		axs[i].plot(idx, x[idx .+ 1]; lw = 0.5, marker = "o", color = "C2")
	end

	axs[i].set_yticks([])
	# axs[i].text(0.99, 0.90, "\$r=$(rs[i])\$"; bbox = bbox,
	# transform = axs[i].transAxes, va = :top, ha=:right, fontsize = 22)
end
# axs[1].set_xlim(0, T)
axs[1].set_xlim(0, T)
add_identifiers!(fig; yloc = 0.8)
fig.tight_layout(pad = 0.25)
fig.subplots_adjust(hspace = 0.1, bottom = 0.1)
# wsave(plotsdir("intermittency"), fig)


# %% Intermittency in logistic map, version with cobweb
using DynamicalSystems, PyPlot
lo = Systems.logistic(0.4)

r3 = 1 + sqrt(8)
rs = [r3, r3 - 0.0001]
# rs = [r3, r3 - 0.0001, r3 - 0.001]
# rs = range(r3; step = -0.0001, length = 3)


chaoticphases = (
[],
[27:42, 141:152, 170:200],
# [0:6, 33:50, 55:76, 83:89, 113:145, 161:200],
)

fig = figure()
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 3; sharex = ax1)
axs = (ax1, ax2)
ax3 = fig.add_subplot(1, 2, 2)
add_identifiers!(fig)

# add inset
ax4 = ax3.inset_axes([0.3, 0.05, 0.4, 0.4])
# ax4.axis("off")

T = 200
dt = 0.05
for i in 1:length(rs)
	set_parameter!(lo, 1, rs[i])
	x = trajectory(lo, T; Ttr = 100)
	axs[i].plot(x; lw = 0.5, marker = "o")
	for idx in chaoticphases[i]
		axs[i].plot(idx, x[idx .+ 1]; lw = 0.5, marker = "o", color = "C2")
	end

	axs[i].set_yticks([])
end
# axs[1].set_xlim(0, T)
axs[1].set_xlim(0, T)

# Do the cobweb plot
function cobweb(t) # transform timeseries t into cobweb (point2D)
    cx = Float64[]; cy = Float64[]
    for i ∈ 1:length(t)-4
        push!(cx, t[i]); push!(cy, t[i])
		push!(cx, t[i]); push!(cy, t[i+1])
    end
    return cx, cy
end

for ax in (ax3, ax4)
	ax.plot([0, 1], [0,1]; lw = 1.0, color = "k")
	xs = 0:0.001:1
	f3 = lo.f.(xs, [rs[end]], 0)
	ax.plot(xs, f3; label = "f")
	f3 = lo.f.(f3, [rs[end]], 0)
	f3 = lo.f.(f3, [rs[end]], 0)
	ax.plot(xs, f3; label = "f³")
	ax.set_xlim(0, 1)
	ax.set_ylim(0, 1)
	x = trajectory(lo, T; Ttr = 100)
	cx, cy = cobweb(x)
	ax.plot(cx, cy; lw = 1.0, alpha = 0.5, color = COLORS[1])
end

# zbox = ((x1, y1), (x2, y2))

fig.tight_layout(pad = 0.25)
fig.subplots_adjust(hspace = 0.1, bottom = 0.1)
# wsave(plotsdir("intermittency"), fig)

# %% Intermittency in Pomeau Maneville map
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

rs = [0.7, 1.0]
Z = 2.6
lo = DiscreteDynamicalSystem(pm_f, 0.3, [rs[1], Z])

fig = figure()
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 3; sharex = ax1)
axs = (ax1, ax2)
ax3 = fig.add_subplot(1, 2, 2)
add_identifiers!(fig)

T = 200
xs = -1:0.001:1

chaoticphases = (
[],
[1:3, 60:83, 160:175],
)


for i in 1:length(rs)
	set_parameter!(lo, 1, rs[i])
	x = trajectory(lo, T; Ttr = 100)
	axs[i].plot(x; lw = 0.5, marker = "o")
	for idx in chaoticphases[i]
		axs[i].plot(idx, x[idx .+ 1]; lw = 0.5, marker = "o", color = "C2")
	end
	axs[i].set_ylim(-1, 1)
	axs[i].set_yticks([])
	# axs[i].text(0.99, 0.90, "\$r=$(rs[i])\$"; bbox = bbox,
	# transform = axs[i].transAxes, va = :top, ha=:right, fontsize = 22)
end
# axs[1].set_xlim(0, T)
axs[1].set_xlim(0, T)

# Do the cobweb plot
function cobweb(t) # transform timeseries t into cobweb (point2D)
    cx = Float64[]; cy = Float64[]
    for i ∈ 1:length(t)-1
        push!(cx, t[i]); push!(cy, t[i])
		push!(cx, t[i]); push!(cy, t[i+1])
    end
    return cx, cy
end


ax3.set_xlim(-0.5, 1.5)
ax4 = ax3.inset_axes([0.4, 0.05, 0.5, 0.5])

for ax in (ax3, ax4)
	f = lo.f.(xs, Ref([rs[end], Z]), 0)
	ax.plot(xs, f, color = "C1")
	ax.plot([-1, 1], [-1,1]; lw = 1.0, color = "k")
	x = trajectory(lo, T; Ttr = 100)
	cx, cy = cobweb(x)
	ax.plot(cx, cy; lw = 1.0, alpha = 0.5, color = COLORS[1])
end

z = 0.04
w = 0.05
zbox = ((z, z), (z+w, z+w))
axis_zoomin!(ax4, ax3, zbox, zbox, "C3")
ax4.set_ylim(z, z+w)
ax4.set_xlim(z, z+w)
ax4.set_xticks([])
ax4.set_yticks([])
ax3.grid(false)

fig.tight_layout(pad = 0.3)
fig.subplots_adjust(hspace = 0.1, bottom = 0.1)
# wsave(plotsdir("intermittency_pm"), fig)

# %% Intermittency in standard map
using DynamicalSystems, PyPlot
ds = Systems.standardmap(;k = 1.2)
u = [2.2, 2.4]

xt = 0:10000:50000
xtl = ["10$(n)" for n in ["⁰", '¹', '²', '³', '⁴', '⁵']]

tr = trajectory(ds, 50000, u)
fig, axs = subplots(1, 2; sharey = true)
θ, p = columns(tr)
n = 0:length(p)-1

axs[1].plot(n[1:1:end], p[1:1:end]; marker = "o", ms = 2, lw = 0.0)
scp = axs[2].scatter(θ, p;
	s = 2, c = n, alpha = 1, edgecolors = "0.5", linewidths = 0.0, cmap = :binary_r,
)

add_identifiers!(fig)

cb = colorbar(scp; ticks = xt)
cb.ax.set_yticklabels(xtl)
cb.set_label("\$n\$", labelpad = 15, rotation = 0)

# intermittency ranges explicit plot
r1 = 1:2700
r2 = 32080:33310
cs = ["C1", "C2"]
for (i, r) in enumerate((r1, r2))
	axs[1].plot(n[r], p[r]; marker = "o", ms = 2, lw = 0.5, color = cs[i])
	axs[2].scatter(θ[r], p[r]; marker = "o", s = 6, color = cs[i])
end
axs[1].set_xticks(xt)
axs[1].set_xticklabels(xtl)
axs[1].set_ylabel("\$p\$")
axs[1].set_xlabel("\$n\$", labelpad = -15)
axs[2].set_xlabel("\$\\theta\$", labelpad = -15)

fig.tight_layout(pad = 0.25)

matplotlib.rcParams["agg.path.chunksize"] = 10000

wsave(plotsdir("intermittency_sm"), fig)

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
wsave(plotsdir("intermittency_manneville"), fig)
