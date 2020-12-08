function naivehist(x::AbstractVector, ε::Real)
    mi, ma = extrema(x)
    hs = zeros(Int((ma-mi)÷ε) + 1)
    is = zeros(Int, length(x))
    for (j, e) in enumerate(x)
        xi = Int((e-mi)÷ε + 1)
        hs[xi] += 1
        is[j] = xi
    end
    return hs ./ length(x), is
end
function naivehist(x::AbstractVector, y::AbstractVector, ε::Real)
    @assert length(x) == length(y)
    x₋, x₊ = extrema(x)
    y₋, y₊ = extrema(y)
    pxy = zeros(Int((x₊ - x₋)÷ε) + 1, Int((y₊ - y₋)÷ε) + 1)
    return naivehist!(pxy, x, y, ε)
end
function naivehist!(pxy, x, y, ε)
    pxy .= 0.0
    x₋ = minimum(x)
    y₋ = minimum(y)
    for j in 1:length(x)
        xv, yv = x[j], y[j]
        xi = Int((xv - x₋)÷ε + 1)
        yi = Int((yv - y₋)÷ε + 1)
        pxy[xi, yi] += 1
    end
    return pxy ./= length(x)
end

function mutualinfo(x, y, ε; base = MathConstants.e)
    px, ix = naivehist(x, ε)
    py, iy = naivehist(y, ε)
    pxy = naivehist(x, y, ε)
    mutualinfo(x, y, px, ix, py, iy, pxy)
end
function mutualinfo(x, y, px, ix, py, iy, pxy, base)
    mi = 0.0
    x₋, x₊ = extrema(x)
    y₋, y₊ = extrema(y)
    for (m, xv) in enumerate(x)
        xi = ix[m]
        for (n, yv) in enumerate(y)
            yi = iy[n]
            ζ = pxy[xi, yi]
            ζ == 0 && continue
            mi += ζ * log(base, ζ/(px[xi]*py[yi]))
        end
    end
    return mi
end

using Random: shuffle!
function mutualinfoshuffle(x, y, ε; M = 10_000, base = MathConstants.e)
    px, ix = naivehist(x, ε)
    py, iy = naivehist(y, ε)
    pxy = naivehist(x, y, ε)
    return mutualinfoshuffle!(copy(x), copy(y), px, ix, py, iy, pxy, M, base)
end
function mutualinfoshuffle!(x, y, px, ix, py, iy, pxy, M, base)
    mi = mutualinfo(x, y, px, ix, py, iy, pxy)
    mis = zeros(M)
    for k in 1:M
        rx = randperm(length(x))
        x = x[rx]; ix = ix[rx]
        ry = randperm(length(x))
        y = y[ry]; iy = iy[ry]
        naivehist!(pxy, x, y, ε)
        mis[k] = mutualinfo(x, y, px, ix, py, iy, pxy, base)
    end
    return mi, mis
end
