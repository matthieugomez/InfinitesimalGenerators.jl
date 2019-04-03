#========================================================================================

Compute the foward looking operator

v_0 * g - ∂(v1 * g) + 0.5 * ∂∂(v2 * g)

========================================================================================#

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

function build_operator(x, v0, v1, v2)
    Δ = make_Δ(x)
    n = length(x)
    T = BandedMatrix(Zeros(n, n), (1, 1))
    build_operator!(T, Δ, v0, v1, v2)
end

function build_operator!(T, Δ, v0, v1, v2)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(T, 0.0)
    # construct matrix T. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    for i in 1:n
        if v1[i] >= 0
            T[min(i + 1, n), i] += v1[i] * invΔxp[i]
            T[i, i] -= v1[i] * invΔxp[i]
        else
            T[i, i] += v1[i] * invΔxm[i]
            T[max(i - 1, 1), i] -= v1[i] * invΔxm[i]
        end
        T[max(i - 1, 1), i] += v2[i] * invΔxm[i] * invΔx[i]
        T[i, i] -= v2[i] * 2 * invΔxm[i] * invΔxp[i]
        T[min(i + 1, n), i] += v2[i] * invΔxp[i] * invΔx[i]
        # Make sure each column sums to zero. Important in some cases: for isntance, otherwise cannot find sdf decomposition in GP model
        T[i, i] += v0[i] - sum(view(T, :, i))
    end
    return T
end


