using DelimitedFiles
cd(@__DIR__)

# 1 = Mackie Glass
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

# 2 = Roessler (1st parameter set)
using DynamicalSystemsBase
roe = Systems.roessler([1.0, 0, 0]; a=0.2, b=0.2, c=5.7)
sroe = trajectory(roe, 500; dt = 0.05, Ttr = 100.0)
writedlm("2.csv", sroe)

# 3 = Lorenz (1st parameter set, dense time sampling)
lo = Systems.lorenz([0, 10.0, 0.0]; σ=10, ρ=28, β=8/3)
tr = trajectory(lo, 100; dt = 0.01, Ttr = 100)
writedlm("3.csv", tr)

# 4 = Roessler, low sampling, 2nd parameter set
ro = Systems.roessler([1.0, 1.0, 1.0], a=0.15, b = 0.2, c=10)
tr = trajectory(ro, 1250; dt = 0.2, Ttr = 100)
writedlm("4.csv", tr)
