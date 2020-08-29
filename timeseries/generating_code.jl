using DelimitedFiles
cd(@__DIR__)

# %% 1 = Mackie Glass timeseries
using DelayDiffEq
function mackey_glass(du,u,h,p,t)
    beta,n,gamma,tau = p
    hist = h(p, t-tau)[1]
    du[1] = (beta*hist)/(1+hist^n) - gamma * u[1]
end
h(p,t) = 0
tau_d = 44
n = 10
β = 0.2
γ = 0.1
δt = 0.5
p = (β,n,γ,tau_d)

tspan = (0.0, 12000.0)
u0 = [1.0]

prob = DDEProblem(mackey_glass,u0,h,tspan,p; constant_lags=tau_d)
alg = MethodOfSteps(Tsit5())
sol = solve(prob,alg; adaptive=false, dt=δt)

s = [u[1] for u in sol.u]
s = s[4001:end]
writedlm("1.csv", s)

# %% 2 = Roessler trajectory (1st parameter set)
using DynamicalSystemsBase
roe = Systems.roessler([1.0, 0, 0]; a=0.2, b=0.2, c=5.7)
sroe = trajectory(roe, 500; dt = 0.05, Ttr = 100.0)
writedlm("2.csv", sroe)

# %% 3 = Lorenz trajectory (1st parameter set, dense time sampling)
lo = Systems.lorenz([0, 10.0, 0.0]; σ=10, ρ=28, β=8/3)
tr = trajectory(lo, 100; dt = 0.01, Ttr = 100)
writedlm("3.csv", tr)

# %% 4 = Roessler trajectory, low sampling, 2nd parameter set
ro = Systems.roessler([1.0, 1.0, 1.0], a=0.15, b = 0.2, c=10)
tr = trajectory(ro, 1250; dt = 0.2, Ttr = 100)
writedlm("4.csv", tr)

# %% 5 = collection of periodic and quasiperiodic orbits
# First part from Henon Heiles
using DynamicalSystems, OrdinaryDiffEq, PyPlot, DelimitedFiles
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
hh = Systems.henonheiles()

xs = Vector{Float64}[]
for u in (
        [0.0, 0.1, 0.5, 0.0], # quasiperiodic
        [0.0, -0.0910355, 0.459522, -0.173339], # higher order quasiperiodic
        [0.0, 0.30266568064689636, 0.42056546332015626, 0.0] # periodic
    )

    tr = trajectory(hh, 500.0, u; Ttr = 10, dt = 0.1, diffeq...)
    x = regularize(tr[:, 1])
    x .+= 0.05randn(length(x))
    push!(xs, x)
end

# Also add a quasi-periodic orbit from Hose Noover for good measure
nh = Systems.nosehoover([0, 1.0, 0])
tr = trajectory(nh, 500.0; Ttr = 10, dt = 0.1, diffeq...)
x = regularize(tr[:, 1])
x .+= 0.05randn(length(x))
push!(xs, x)

# Second part is chaotic & period-4 roessler
ro = Systems.roessler()
for c in (4.0, 5.7)
    set_parameter!(ro, 3, c)
    tr = trajectory(ro, 500.0; Ttr = 100, dt = 0.1, diffeq...)
    x = regularize(tr[:, 1])
    x .+= 0.05randn(length(x))
    push!(xs, x)
end

# Third part from standard map
sm = Systems.standardmap()

for (i, u) in enumerate((
        SVector(0.8121, 1.6243),
        SVector(0.877, 1.565),
    ))
    tr = trajectory(sm, 5000, u; Ttr = 10)
    x = regularize(tr[:, 1])
    i == 1 && (x .+= 0.05randn(length(x)))
    push!(xs, x)
end

# Fourth part is just oscillatory ARMA, a-la chapter 7
using Random, Statistics
Random.seed!(77163)
η = randn(15000)
s = ones(15000)
for n in 5:15000
    s[n] = 1.625s[n-1] - 0.284s[n-2] - 0.355s[n-3] + η[n] - 0.26η[n-1]
end
x = s[end-5000:end]
x ./= std(x)
push!(xs, x)

using Statistics, FFTW
# for x in xs
#     plot(x[1:1000], alpha = 0.5)
#     r = x
#     P = abs2.(rfft(r .- mean(r)))
#     f = FFTW.rfftfreq(length(r))
#     plot(f, P)
# end

writedlm("5.csv", hcat(xs...))

# %% 6 = towel map trajectory
using DynamicalSystems, DelimitedFiles

ds = Systems.towel()

tr = trajectory(ds, 10000; Ttr = 100)

writedlm("6.csv", tr)
