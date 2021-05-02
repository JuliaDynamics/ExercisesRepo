using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("colorscheme.jl"))
using DynamicalSystems, InteractiveDynamics
import GLMakie
using AbstractPlotting
using LinearAlgebra
using Statistics

# TODO: Add rotation

# Code for 3D animation
ds = Systems.lorenz()
u0 = trajectory(ds, 10; Ttr = 100)[end]
dt = 0.005
N = 100
framerate = 30
systemtitle = "Lorenz63"

ds = Systems.towel()
u0 = trajectory(ds, 10; Ttr = 100)[end]
dt = 1
N = 15 # total steps to take
framerate = 2
systemtitle = "towelmap"

# %%
tinteg = tangent_integrator(ds, 3; u0)
Y = get_deviations(tinteg)

# Now I define functions that make an ellipsoid given the point vectors
function ellipsoid_transform(Y)
    svddecomp = svd(Y)
    d = svddecomp.S
    s = d ./ maximum(d)
    Q = svddecomp.U * svddecomp.Vt # rotation matrix
    θx = atan(Q[3,2],Q[3,3])
    θy = atan(-Q[3,1], sqrt(Q[3,2]^2+Q[3,3]^2))
    θz = atan(Q[2,1], Q[1,1])
    M = AbstractPlotting.rotationmatrix_x(θx) *
        AbstractPlotting.rotationmatrix_y(θy) *
        AbstractPlotting.rotationmatrix_z(θz) *
        AbstractPlotting.scalematrix(Vec3f0(s...))
    return M, s
end
function transform_mesh(msh, mat4x4)
    pos_trans = Point3f0.(Ref(mat4x4) .* to_ndim.(Point4f0, spheremesh.position, 0))
    AbstractPlotting.GeometryBasics.Mesh(pos_trans, AbstractPlotting.GeometryBasics.faces(msh))
end
function ellipsoid_from_pointcloud(v, msh = spheremesh)
    M, s = ellipsoid_transform(v)
    return transform_mesh(msh, M), s
end
spheremesh = AbstractPlotting.GeometryBasics.mesh(Sphere(Point3f0(0), 1))

ellipsoid, s = ellipsoid_from_pointcloud(Y, spheremesh)
colors = [RGBAf0((abs.(u)./norm(u))..., 1) for u in spheremesh.position]

ellipsobs = Observable(ellipsoid)

# Setup figure and plot everything

# fig = Figure(resolution = (1200, 600))
fig = Figure(resolution = (1200, 600)); display(fig)
ax3D = Axis3(fig[1, 1]; aspect = :data)
mesh!(ellipsobs, color = colors)

ax3D.xticklabelsvisible = true
ax3D.yticklabelsvisible = true
ax3D.zticklabelsvisible = true
ax3D.azimuth = 3.92
ax3D.elevation = 0.38
ax3D.title = rpad(systemtitle *", t = 0", 40, ' ')

ax2 = Axis(fig[1, 2]; width = 400) # this axis has values of principal lengths
ax2.title = "ellipsoid axes"
sobs = Observable(s)
barplot!(ax2, 1:3, sobs; color = to_color.(COLORSCHEME[1:3]))

display(fig)

# %%
t = 0
record(
        fig, joinpath(@__DIR__, "volume_growth_$(systemtitle).mp4"), 1:N;
        framerate
        ) do i
    DynamicalSystems.step!(tinteg, dt)
    Y = get_deviations(tinteg)
    ellipsobs[], sobs[] = ellipsoid_from_pointcloud(Y)
    # s[] = principals(newus)
    # autolimits!(ax2)
    # must set new limits centered on u
    # maxε = maximum(norm(u) for u in newus)
    # maxε = maximum(s)
    # ax3D.limits = map(j -> (-maxε, +maxε), (1,2,3))
    global t = round(tinteg.t; digits = 5)
    ax3D.title = rpad(systemtitle *", t = $(t)", 40, ' ')
    ax3D.elevation = 0.38
    ax3D.azimuth = 3.92
end
