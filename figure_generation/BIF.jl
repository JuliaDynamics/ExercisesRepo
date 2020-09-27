include("style.jl")

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

fig = figure(figsize = (figx, 1.1figx))
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
# fsave(joinpath(figdir, "orbit_diagrams"), fig; pdf = false, dpi = 600)

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
# iter = BK.PALCIterable(F, Jac_m, [0.8], 1., (BK.@lens _), opts; verbosity = 2)
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
xlabel("\$\\epsilon\$")
ylabel("\$T^*\$")
tight_layout()

fig.subplots_adjust(left = 0.18, bottom = 0.2, top = 0.97)
# fsave(joinpath(figdir, "bif_example"), fig)
