
struct Derivative
    x::Vector{Float64}
end

function FirstDerivative(x, y::T, i;  bc = (0, 0), direction = :up)
    if method == :up
        Δxp = x[min(i, n-1)+1] - x[min(i, n-1)]
        (i < length(x)) ? (y[i+1] - y[i]) / Δxp : convert(T, bc[end])
    elseif method == :down
        Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
        (i > 1) ? (y[i] - y[i-1]) / Δxm : convert(T, bc[1])
    end
end

function SecondDerivative(x, y::T, i, bc = (0, 0))
    Δxp = x[min(i, n-1)+1] - x[min(i, n-1)]
    Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
    Δx = (Δxm + Δxp) / 2
    (1 < i < length(x)) ? (y[i + 1] / (Δxp * Δx) + y[i - 1] / (Δxm * Δx) - 2 * y[i] / (Δxp * Δxm)) : ((i == 1) ? (y[2] / (Δxp * Δx) + (y[1] - bc[1] * Δxm) / (Δxm * Δx) - 2 * y[1] / (Δxp * Δxm)) : ((y[end] + bc[end] * Δxp) / (Δxp * Δx) + y[end - 1] / (Δxm * Δx) - 2 * y[end] / (Δxp * Δxm)))


    if method == :upwind
        Δxp = x[min(i, n-1)+1] - x[min(i, n-1)]
        Δfp = y[min(i, n-1)+1] - y[min(i, n-1)]
        if Δxp != 0
            return Δfxp / Δxp
        else
            return 0.0
        end
    else
        Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
        Δfm = f[max(i-1, 1) + 1] - f[max(i-1, 1)]
        if Δxm != 0
            return Δfxm / Δxm
        else
            return 0.0
        end
    end
end
    elseif method == :down
        return x[min(i, n-1)+1] - x[min(i, n-1)]
struct ∂down
    grid::Vector{Float64}
end

struct ∂down
    grid::Vector{Float64}
end
upwinding(wgrid)