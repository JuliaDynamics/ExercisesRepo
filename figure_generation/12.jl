using DrWatson
@quickactivate "ExercisesRepo"
include(srcdir("style.jl"))
using DynamicalSystems, PyPlot, LinearAlgebra, Statistics

# %% Sparse state space
function sparse_f!(xnew, x, p, n)
    xnew[1] = 4*x[1]*(1-x[1])
    @inbounds for i in 2:length(x)
        xnew[i] = p[i]*x[i-1]
    end
    return
end

function sparse_jacob!(J, x, p, n)
    J[1,1] = 4 - 8x[1]
end

function make_sparse_system(D, p = [0.2*rand() + 0.87 for _ in 1:D-1])
    u0 = rand(D)
    # Optimized Jacobian in the form of tri-diagonal
    J = Tridiagonal(zeros(D-1), zeros(D), zeros(D-1))
    # Add entries of linear components, that will not change in time
    for i in 2:D; J[i, i-1] = p[i-1]; end
    # Adjust the first entry of the Jacobian, which depends on the state
    sparse_jacob!(J, u0, p, 0)
    sparse = DiscreteDynamicalSystem(sparse_f!, u0, p, sparse_jacob!, J)
end

sparse = make_sparse_system(400)
# A = trajectory(sparse, 1000; save_idxs = 1)
#
# plot(A[:, 1])

# Ds = 10:20:400
#
# λs = lyapunovspectrum(sparse, 10000, 10)
#
# Δs[i] = kaplanyorke_dim(sort!(λs; rev = true))

# data 1: kaplan yorke for increasing D
fig, axs = subplots(1,2; figsize = (figx, figy))

Ds = 10:10:200
Δs = zeros(length(Ds))

# TODO: add produce_or_load here
for (i, D) in enumerate(Ds)
    @show D
    sparse = make_sparse_system(D)
    @time λs = lyapunovspectrum(sparse, 10000, 10)
    Δs[i] = kaplanyorke_dim(sort!(λs; rev = true))
    break
end
# plot 1


axs[1].plot(Ds, Δs)
axs[1].set_xlabel("\$D\$"; labelpad = -10)
axs[1].set_ylabel("\$\\Delta^{(L)}\$"; labelpad = -5)
axs[1].set_yticks(1:4)
axs[1].set_xticks(10:40:200)


# %% data 2: perturbation growth lala
# fig, axs = subplots(1,2; figsize = (figx, figy))
using Random
# Random.seed!(77670011)
Ds = 10 .^ (1:4)
for (d, D) in enumerate(Ds)
    @show D
    sparse = make_sparse_system(D)


    S = 1000 # samples
    n = 50 # max iteration
    lD = zeros(n, S)

    for j in 1:S
        Q0 = normalize!(rand(D, 1))
        while any(isnan, Q0)
            Q0 = normalize!(rand(D, 1))
        end
        Q0 .* 1e-6
        tinteg = tangent_integrator(sparse, Q0; u0 = rand(D))
        g = zeros(n)
        for i in 1:n
            step!(tinteg)
            g[i] = norm(get_deviations(tinteg))
        end

        # create instantaneous exponential growth rates
        lD[:, j] = log.(g) ./ (1:n)
    end

    lDμ = vec(mean(lD; dims = 2))
    lDσ = vec(std(lD; dims = 2))
    axs[2].plot(1:n, lDμ; color = "C$(d)", label = "\$D=$D\$")
    # axs[2].fill_between(1:n, lDμ .- lDσ, lDμ .+ lDσ; color = "C$(d)", alpha = 0.5)
end

# add lyapunov exponent indication
axs[2].axhline(log(2); ls = "dashed", color = "C0")
axs[2].set_xlabel("\$n\$"; labelpad = -10)
axs[2].set_xticks(0:16:48)
axs[2].set_ylabel("\$\\lambda_\\mathrm{local}\$")
axs[2].legend()
axs[2].set_ylim(-0.1, 0.75)
add_identifiers!(fig)
fig.tight_layout(;pad = 0.5)
wsave(plotsdir("sparse_space"), fig)


# %% Rate-dependent example
using DynamicalSystems, PyPlot, Roots, ForwardDiff
fig = figure(figsize = (figx/2, figy))
axbif = subplot(1,1,1)

α(T) = 0.5 - 0.2*tanh((T-263)/4)
dTdt(T, ε=0.65) = 1 - α(T) - 10*ε*(0.002T)^4
es = 0.3:0.005:0.9
for e in es
    f = (T) -> dTdt(T, e)
    roots = Roots.find_zeros(f, 50, 500)
    for (i, r) in enumerate(roots)
        c = ForwardDiff.derivative(f, r) < 0 ? "k" : "w"
        axbif.plot(e, r, marker = "o", mec = "k", mew = 1,
        markersize = 6, mfc = c )
    end
end


axe = axbif.inset_axes(;bounds = [0.6, 0.6, 0.35, 0.35], transform = axbif.transAxes)
axbif.set_ylim(210, 370)

function rate_dependent_system(du, u, p, t)
    Δε, Δt, t0 = p
    T, ε = u
    dεdt = (-(t - t0)*Δε / Δt^2) * exp(-(t-t0)^2/(2*Δt^2))
    du[1] = dTdt(T, ε)
    du[2] = dεdt
    nothing
end

t0 = 2500.0

for (i, p) in enumerate([(0.1, 40.0), (0.15, 40.0), (0.1, 60.0)])
    p0 = [p..., t0] # if Δε is 0, it doesn't matter
    u0 = [278, 0.78] # lol gotta rework this because ε goes beyond 1.

    ds = ContinuousDynamicalSystem(rate_dependent_system, u0, p0)

    maxt = 5000
    tr = trajectory(ds, maxt; dt = 1.0)
    T, ε = columns(tr)
    axbif.plot(ε, T; color = "C$(i)")
    axe.plot(0:maxt, ε; color = "C$(i)")
end
axe.set_xlim(2300, 2700)
axe.set_xlabel("\$t\$")
axe.set_xticks([])
axe.set_ylabel("\$\\epsilon\$"; labelpad = -20)
axe.grid(false)
axbif.set_xticks(0.3:0.2:0.9)
axbif.set_yticks(210:50:370)
axbif.set_xlabel("\$\\epsilon\$"; labelpad = -20)
axbif.set_ylabel("\$T\$"; labelpad = -20)
fig.tight_layout(pad = 0.25)
wsave(plotsdir("rate_dependent"), fig)

# %% Basin stability of magnetic pendulum
using DynamicalSystems, PyPlot, LinearAlgebra, Statistics

struct MagneticPendulum{T<:AbstractFloat}
    magnets::Vector{SVector{2, T}}
end

function (m::MagneticPendulum)(u, p, t)
    x, y, vx, vy = u
    γ1, γ2, γ3, d, α, ω = p
    γ = (γ1, γ2, γ3)
    dx, dy = vx, vy
    dvx, dvy = @. -ω^2*(x, y) - α*(vx, vy)
    for (i, ma) in enumerate(m.magnets)
        δx, δy = (x - ma[1]), (y - ma[2])
        D = sqrt(δx^2 + δy^2 + d^2)
        dvx -= γ[i]*(x - ma[1])/D^3
        dvy -= γ[i]*(y - ma[2])/D^3
    end
    return SVector(dx, dy, dvx, dvy)
end

p0 = [1, 1, 0.2, 0.3, 0.2, 0.5]
u0 = [sincos(0.12553*2π)..., 0, 0]
m = MagneticPendulum([SVector(cos(2π*i/3), sin(2π*i/3)) for i in 1:3])
ma = ContinuousDynamicalSystem(m, u0, p0)

using LinearAlgebra

const xmin, xmax = -3, 3

function high_quality_statespace(ma, g = range(xmin, xmax; length = 200))
    ∫ = integrator(ma)

    c = fill(0, length(g), length(g))
    t = similar(c, Float64)

    for (i, x) ∈ enumerate(g)
        # @show x
        for (j, y) in enumerate(g)
            reinit!(∫, SVector(x, y, 0, 0))
            t0 = ∫.t0
            step!(∫, 100.0)
            while ∫.u[3]^2 + ∫.u[4]^2 > 1e-3
                step!(∫)
            end
            s = SVector(∫.u[1], ∫.u[2])
            dmin, k = findmin([(s-m)⋅(s-m) for m in ma.f.magnets])
            # t[i,j] = ∫.t - ∫.t0
            c[i,j] = k
            # break
        end
    end
    return c
end

function random_statespace_fraction(∫, T = 10000)
    c = 0 # count of states converged to attractor of interest
    for i in 1:T
        x = rand()*(xmax-xmin) + xmin
        y = rand()*(xmax-xmin) + xmin
        reinit!(∫, SVector(x, y, 0, 0))
        step!(∫, 100.0)
        while ∫.u[3]^2 + ∫.u[4]^2 > 1e-3
            step!(∫)
        end
        s = SVector(∫.u[1], ∫.u[2])
        dmin, k = findmin([(s-m)⋅(s-m) for m in ma.f.magnets])
        k == 3 && (c += 1)
    end
    return c/T
end


γs = 0:0.01:1
Fs = zero(γs)
Fσs = zero(γs)
λs = zero(γs)
∫ = integrator(ma)

for (i, γ) in enumerate(γs)
    @show (i, γ)
    p0[3] = γ
    allfracs = [random_statespace_fraction(∫, 1000) for i in 1:10]
    Fs[i] = mean(allfracs)
    Fσs[i] = std(allfracs)
    J = Array(ma.jacobian(SVector(m.magnets[3]..., 0, 0), p0, 0))
    eee = eigvals(J)
    λs[i] = maximum(real.(eee))
end


# %%
fig, axs = subplots(1,3; figsize = (figx, figy))
axs[1].plot(γs, Fs; label = "\$F\$")
axs[1].fill_between(γs, Fs .- Fσs, Fs .+ Fσs; color = "C0", alpha = 0.5)
# axs[1].errorbar(γs, Fs; label = "\$F\$", yerr = Fσs)
# axs[1].plot(γs, Fs; label = "\$F\$")
axs[1].plot(γs, λs; label = "\$\\lambda_1\$")
axs[1].legend()
axs[1].set_xlabel("\$\\gamma_3\$"; labelpad = -20)
axs[1].set_xticks([0, 1])
axs[1].set_xlim(0,1)
axs[1].set_yticks(-0.1:0.1:0.3)

# perform 2 nice detailed plots
LC =  matplotlib.colors.ListedColormap
cmap = LC([matplotlib.colors.to_rgb("C$k") for k in 0:2])
γplots = (0.7, 0.25)
for (i, γ) in enumerate(γplots)
    p0[3] = γ
    g = range(xmin, xmax; length = 1000)
    c = high_quality_statespace(ma, g)
    axs[i+1].pcolormesh(g, g, c'; cmap = cmap, shading = "gouraud")
    # axs[i+1].pcolormesh(g, g, c'; cmap = cmap)
end

for i in 2:3;
    axs[i].set_xlabel("\$x\$")
    axs[i].set_ylabel("\$y\$")
    axs[i].set_xticks([])
    axs[i].set_yticks([])
end

add_identifiers!(fig)
fig.tight_layout(;pad = 0.5)
wsave(plotsdir("basin_stability"), fig)

# %% phase induced transitions
