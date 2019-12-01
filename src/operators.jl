
#========================================================================================

Compute the operator
𝔸f = v_0 * f + v1 * ∂(f) + v2 * ∂∂(f)
Note that
𝔸'g = v_0 * g - ∂(v1 * g) + ∂∂(v2 * g)

========================================================================================#

function operator(x::AbstractVector, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    n = length(x)
    𝔸 = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    Δ = make_Δ(x)
    operator!(𝔸, Δ, v0, v1, v2)
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

function operator!(𝔸, Δ, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    # The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(𝔸, 0)
    @inbounds for i in 1:n
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
    end
    c = sum(𝔸, dims = 1)
    for i in 1:n
        𝔸[i, i] += v0[i] - c[i]
    end
    return adjoint(𝔸)
end



#========================================================================================

Compute the principal eigenvector and eigenvalue of 𝔸
By definition, it is the one associated with a positive eigenvector.
In particular, it must be real.

B = -𝔸 is a Z matrix (all off diagonal are negative). Therefore, there exists a positive s such that sI + A has all positive entries. Applying Perron Frobenus, there a unique largest eigenvalue for sI + A, which is real, and the correspongind eigenctor is strictly positive.
Note that, in particular, it is the eigenvalue with largest real part, which means that I can look for the eigenvalue with largest real part 



If, moreover, B, is a M-matrix, then all its eigenvalues have positive real part. Therefore, all the eigenvalues of A have negative real part. Therefore, the eigenvalue with largest real part is also the eigenvalue with smallest magnitude.

========================================================================================#
function principal_eigenvalue(𝔸::AbstractMatrix; which = :SM, eigenvector = :right)
    f, η, g = nothing, nothing, nothing
    if which == :SM
        if eigenvector ∈ (:right, :both)
            vals, vecs = Arpack.eigs(𝔸; nev = 1, which = :SM)
                η = vals[1]
                f = vecs[:, 1]
        end
        if eigenvector ∈ (:left, :both)
            vals, vecs = Arpack.eigs(adjoint(𝔸); nev = 1, which = :SM)
            η = vals[1]
            g = vecs[:, 1]
        end 
    elseif which == :LR
        # While Arpack accepts SM, it often fails. Moreover it does not give the "right" eigenvector in term of multiplicity.
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
    end
    return clean_eigenvector_left(g), clean_eigenvalue(η), clean_eigenvector_right(f)
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
    g = abs.(g)
    g ./ sum(g)
end

clean_eigenvector_right(::Nothing) = nothing
clean_eigenvector_right(f::Vector) = abs.(f) / sum(abs.(f)) .* length(f)
