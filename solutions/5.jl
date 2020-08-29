# 5.5
using DynamicalSystems, PyPlot, Statistics

x = readdlm("7.csv")
W = 100 # timeseries length
w = 20
@assert w < W

L = 1:w:length(x)-W
os = (3, 4, 5) # orders

H = zeros(length(L), length(os))

for (j, i) in enumerate(L)
    for (k, o) in enumerate(os)
        H[j, k] = permentropy(view(x, i:i+W), o)
    end
end

for (k, o) in enumerate(os)
    H[:, k] ./= std(view(H, :, k))
end

figure(); plot(H)
