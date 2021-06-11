# Magnetic pendulum basins of attraction
using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, Random

α=0.2; ω=1.0; d=0.3
gx = gy = range(-5, 5; length = 1500)
config = @strdict(α, ω, d, gx, gy)

function magnetic_basins(config)
    @unpack α, ω, d, gx, gy = config
    ma = Systems.magnetic_pendulum(; α, ω, d)
    @time basins, attractors = basins_general(gx, gy, ma; idxs = 1:2)
    return @strdict(gx, gy, basins, attractors)
end

# produce high quality
produce_or_load(datadir(), config, magnetic_basins; prefix = "magnetic_basins")

# produce zoomed version
gx = range(1.80, 1.95; length = 1000)
gy = range(0, 0.12; length = 1000)
config = @strdict(α, ω, d, gx, gy)
produce_or_load(datadir(), config, magnetic_basins; prefix = "magnetic_basins_zoomed")

# Plot this
gx, gy, basins, attractors = @unpack wload(datadir(savename("magnetic_basins", config)))

LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])

fig = figure(figsize=(figx/2, figy))
pcolormesh(gx, gy, basins'; cmap, shading = "gouraud")
gca().set_aspect("equal")
xticks([-5, 5])
yticks([-5, 5])
xlabel("\$x\$", labelpad=-30)
ylabel("\$y\$", labelpad=-30)
for m in attractors
    m1, m2 = columns(m)
    scatter(m1, m2; color = "white", edgecolor = "black", zorder = 99, s = 100)
end
tight_layout(pad=0.3)
subplots_adjust(bottom = 0.1, top = 0.95)
wsave(plotsdir("magneticpendulum"), fig)

# Plot zoomed version
gx, gy, basins, attractors = @unpack wload(datadir(savename("magnetic_basins_zoomed", config)))
LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])

fig = figure(figsize=(figx/2, figx/2))
pcolormesh(gx, gy, basins'; cmap, shading = "gouraud")
gca().set_aspect("equal")

xticks([1.8, 1.95], size = 20)
yticks([0, 0.12], size = 20)

fig.savefig(joinpath(figdir, "magneticpendulum_zoom.png"))