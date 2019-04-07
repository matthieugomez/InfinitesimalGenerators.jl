#========================================================================================

Compute the operator
𝔸f = v_0 * f + v1 * ∂(f) + 0.5 * v2 * ∂∂(f)
𝔸'g = v_0 * g - ∂(v1 * g) + 0.5 * ∂∂(v2 * g)

========================================================================================#

function operator(x, v0, v1, v2)
    𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    operator!(𝔸, make_Δ(x), v0, v1, v2)
end

function operator!(𝔸, Δ, v0, v1, v2)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(𝔸, 0.0)
    # construct matrix T. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    for i in 1:n
        if v1[i] >= 0
            𝔸[i, min(i + 1, n)] += v1[i] * invΔxp[i]
            𝔸[i, i] -= v1[i] * invΔxp[i]
        else
            𝔸[i, i] += v1[i] * invΔxm[i]
            𝔸[i, max(i - 1, 1)] -= v1[i] * invΔxm[i]
        end
        𝔸[i, max(i - 1, 1)] += v2[i] * invΔxm[i] * invΔx[i]
        𝔸[i, i] -= v2[i] * 2 * invΔxm[i] * invΔxp[i]
        𝔸[i, min(i + 1, n)] += v2[i] * invΔxp[i] * invΔx[i]
    end
    # Make sure each row sums to zero. Important in some cases: for isntance, otherwise cannot find sdf decomposition in GP model
    c = sum(𝔸, dims = 2)
    for i in 1:n
        𝔸[i, i] += v0[i] - c[i]
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

#========================================================================================

Compute the principal eigenvector and eigenvalue of an operator

========================================================================================#
clean_eigenvalue(η::Union{Nothing, Real}) = η

function clean_eigenvalue(η::Complex)
    if abs(imag(η) .>= eps())
        @warn "Principal Eigenvalue has some imaginary part $(η)"
    end
    real(η)
end
clean_density(::Nothing) = nothing
clean_density(v::Vector) = abs.(v) ./ sum(abs.(v))

clean_f(v::Vector) = abs.(v)
clean_f(::Nothing) = nothing

function principal_eigenvalue(T; method = :krylov, eigenvector = :right)
    η = nothing
    if method == :krylov
        vl, η, vr = principal_eigenvalue_krylov(T; eigenvector = eigenvector)
        if η == nothing
            @warn "Krylov Methods Failed"
        end
    end
    if η == nothing
        # use SuiteSparse maybe? LU decomposition sometimes?
        vl, η, vr = principal_eigenvalue_BLAS(convert(Matrix{Float64}, T); eigenvector = eigenvector)
    end
    return clean_density(vl), clean_eigenvalue(η), clean_f(vr)
end

# I could also use Arpack.eigs but it seems slower
function principal_eigenvalue_krylov(T; eigenvector = :right)
    vl, η, vr = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        vals, vecs, info = KrylovKit.eigsolve(T, 1, :LR, maxiter = size(T, 1))
        if info.converged > 0
            η = vals[1]
            vr = vecs[1]
        end
    end
    if eigenvector ∈ (:left, :both)
        vals, vecs, info = KrylovKit.eigsolve(T', 1, :LR, maxiter = size(T, 1))
        if info.converged > 0
            η = vals[1]
            vl = vecs[1]
        end
    end 
    return vl, η, vr
end

function principal_eigenvalue_BLAS(T; eigenvector = :right)
    vl, η, vr = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        e = eigen(T)
        _, out = findmax(real.(e.values))
        η = e.values[out]
        vr = e.vectors[:, out]
    end
    if eigenvector ∈ (:left, :both)
        e = eigen(copy(T'))
        _, out = findmax(real.(e.values))
        η = e.values[out]
        vl = e.vectors[:, out]
    end 
    return vl, η, vr
end
