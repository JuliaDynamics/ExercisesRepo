
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

# wsave(plotsdir("intermittency_sm"), fig)


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
