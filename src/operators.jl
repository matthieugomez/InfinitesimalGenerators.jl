
#========================================================================================

Compute the operator
𝔸f = v_0 * f + v1 * ∂(f) + v2 * ∂∂(f)
Note that
𝔸'g = v_0 * g - ∂(v1 * g) + ∂∂(v2 * g)

========================================================================================#

function operator(x::AbstractVector, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    𝔸 = Tridiagonal(zeros(length(x)-1), zeros(length(x)), zeros(length(x)-1))
    operator!(𝔸, make_Δ(x), v0, v1, v2)
end

function operator!(𝔸::AbstractMatrix, Δ, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
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

Compute the principal eigenvector and eigenvalue of 𝔸

========================================================================================#
function principal_eigenvalue(𝔸::AbstractMatrix; method = :krylov, eigenvector = :right)
    η = nothing
    if method == :krylov
        g, η, f = principal_eigenvalue_krylov(𝔸; eigenvector = eigenvector)
        if η == nothing
            @warn "Krylov Methods Failed"
        end
    end
    if η == nothing
        g, η, f = principal_eigenvalue_BLAS(convert(Matrix{Float64}, 𝔸); eigenvector = eigenvector)
    end
    return clean_eigenvector_left(g), clean_eigenvalue(η), clean_eigenvector_right(f)
end

# I could also use Arpack.eigs but it seems slower
function principal_eigenvalue_krylov(𝔸::AbstractMatrix; eigenvector = :right)
    g, η, f = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        vals, vecs, info = KrylovKit.eigsolve(𝔸, 1, :LR, maxiter = size(𝔸, 1))
        if info.converged > 0
            η = vals[1]
            f = vecs[1]
        end
    end
    if eigenvector ∈ (:left, :both)
        vals, vecs, info = KrylovKit.eigsolve(adjoint(𝔸), 1, :LR, maxiter = size(𝔸, 1))
        if info.converged > 0
            η = vals[1]
            g = vecs[1]
        end
    end 
    return g, η, f
end

function principal_eigenvalue_BLAS(𝔸::AbstractMatrix; eigenvector = :right)
    g, η, f = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        e = eigen(𝔸)
        _, out = findmax(real.(e.values))
        η = e.values[out]
        f = e.vectors[:, out]
    end
    if eigenvector ∈ (:left, :both)
        e = eigen(copy(adjoint(𝔸)))
        _, out = findmax(real.(e.values))
        η = e.values[out]
        g = e.vectors[:, out]
    end 
    return g, η, f
end

clean_eigenvalue(η::Union{Nothing, Real}) = η
function clean_eigenvalue(η::Complex)
    if abs(imag(η) .>= eps())
        @warn "Principal Eigenvalue has some imaginary part $(η)"
    end
    real(η)
end

clean_eigenvector_left(::Nothing) = nothing
function clean_eigenvector_left(g::Vector)
    abs.(g) ./ sum(abs.(g))
end

clean_eigenvector_right(::Nothing) = nothing
clean_eigenvector_right(f::Vector) = abs.(f) / sum(abs.(f)) .* length(f)
