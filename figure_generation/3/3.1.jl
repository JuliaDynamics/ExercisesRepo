using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

# sensitive dependence
Random.seed!(5)
close("all")
fig = figure(figsize=(figx/2, figy))
ax = gca()
u0 = [10,10.0,10]
lo = Systems.lorenz([10,10.0,10])

for i in 1:3
    u = u0 .+ i*1e-3
    tr = trajectory(lo, 15, u)
    plot(0:0.01:15, tr[:, 1],  c = COLORS[[2,3,1][i]])
end

xlabel("\$t\$",labelpad = -10)
yticks(-15:10:15)
xticks(0:5:15)
ylabel("\$x\$",labelpad = -20)
tight_layout(pad = 0.25)
subplots_adjust(bottom = 0.15, left = 0.12)
wsave(plotsdir("sensitive"), fig)
