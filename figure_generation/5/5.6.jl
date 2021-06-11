using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

he = Systems.henon()
N, α = 10000, 2
X = trajectory(he, N, Ttr = 10)
ε = 10 .^ (-5.0:0.1:1)
H = genentropy(α, ε, X)
C = correlationsum(X, ε)

fig, axs = subplots(2, 1; figsize = (figx/2, figx/2), sharex = true)

axs[1].plot(log.(ε), -H, COLORS[1])

x = -5:0
D = 1.22

axs[2].plot(log.(ε), log.(C), COLORS[1])

for (ax, x) in zip(axs, (-5:0, -10:0))
    Dl = D .* x
    Dl = Dl .- 2
    ax.plot(x, Dl, c = COLORS[3])
    ax.axvspan(x[1], x[end], color = COLORS[2], alpha = 0.25)
    ax.grid("on")
    ax.text(x[2] + 1.5, Dl[2], "\$D_2\$ = $D", color = COLORS[3], size = 30)

end

axs[2].set_xticks(-10:2:2)
axs[1].tick_params(labelbottom=false)
axs[2].set_xlabel("\$\\log(\\varepsilon)\$")
axs[1].set_ylabel("\$-H_$α\$", labelpad = -20)
axs[2].set_ylabel("\$\\log(C)\$", labelpad = -20)
axs[1].set_yticks(-10:2:0)
axs[1].set_yticklabels(["-10", "", "", "", "", "0"])
axs[2].set_yticks(0:-2:-12)
axs[2].set_yticklabels(["0", "", "", "",  "","", "-12"])
fig.tight_layout()
fig.subplots_adjust(left=0.14, bottom = 0.15, hspace = 0.1)
# wsave(joinpath(figdir, "dimension"), fig)
