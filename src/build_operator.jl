#========================================================================================

Compute the foward looking operator

𝔸g = v_0 * g - ∂(v1 * g) + 0.5 * ∂∂(v2 * g)
𝔸'f = v_0 * f + v1 * ∂(f) + 0.5 * v2 * ∂∂(f)

========================================================================================#

function build_operator(x, v0, v1, v2)
    𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    build_operator!(𝔸, make_Δ(x), v0, v1, v2)
end

function build_operator!(𝔸, Δ, v0, v1, v2)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(𝔸, 0.0)
    # construct matrix T. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    for i in 1:n
        if v1[i] >= 0
            𝔸[min(i + 1, n), i] += v1[i] * invΔxp[i]
            𝔸[i, i] -= v1[i] * invΔxp[i]
        else
            𝔸[i, i] += v1[i] * invΔxm[i]
            𝔸[max(i - 1, 1), i] -= v1[i] * invΔxm[i]
        end
        𝔸[max(i - 1, 1), i] += v2[i] * invΔxm[i] * invΔx[i]
        𝔸[i, i] -= v2[i] * 2 * invΔxm[i] * invΔxp[i]
        𝔸[min(i + 1, n), i] += v2[i] * invΔxp[i] * invΔx[i]
        # Make sure each column sums to zero. Important in some cases: for isntance, otherwise cannot find sdf decomposition in GP model
        𝔸[i, i] += v0[i] - sum(view(𝔸, :, i))
    end
    return 𝔸
end

function make_Δ(x)
    n = length(x)
    Δxm = zero(x)
    Δxm[1] = x[2] - x[1]
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zero(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δxp[end] = x[n] - x[n-1]
    Δx = (Δxm .+ Δxp) / 2
    return x, 1 ./ Δx, 1 ./ Δxm, 1 ./ Δxp
end
