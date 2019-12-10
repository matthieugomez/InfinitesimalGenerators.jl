
#========================================================================================

Compute the operator
𝔸f = v_0 * f + v1 * ∂(f) + v2 * ∂∂(f)
Note that
𝔸'g = v_0 * g - ∂(v1 * g) + ∂∂(v2 * g)

========================================================================================#
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
    adjoint(𝔸)
end

#========================================================================================

Compute the principal eigenvector and eigenvalue of 𝔸
By definition, it is the one associated with a positive eigenvector.
In particular, it must be real.

B = -𝔸 is a Z matrix (all off diagonal are negative). Therefore, there exists a positive s such that sI + A has all positive entries. Applying Perron Frobenus, there a unique largest eigenvalue for sI + A, which is real, and the correspongind eigenctor is strictly positive.
Note that, in particular, it is the eigenvalue with largest real part, which means that I can look for the eigenvalue with largest real part 



If, moreover, B, is a M-matrix, then all its eigenvalues have positive real part. Therefore, all the eigenvalues of A have negative real part. Therefore, the eigenvalue with largest real part is also the eigenvalue with smallest magnitude.

========================================================================================#
function principal_eigenvalue(𝔸::AbstractMatrix; which = :SM, eigenvector = :right, r0 = ones(size(𝔸, 1)))
    f, η, g = nothing, nothing, nothing
    if which == :SM
        if eigenvector ∈ (:right, :both)
            vals, vecs = Arpack.eigs(𝔸; v0 = r0, nev = 1, which = :SM)
                η = vals[1]
                f = vecs[:, 1]
        end
        if eigenvector ∈ (:left, :both)
            vals, vecs = Arpack.eigs(adjoint(𝔸); nev = 1, which = :SM)
            η = vals[1]
            g = vecs[:, 1]
        end 
    elseif which == :LR
        # Arpack LR tends to fail if the LR is close to zero, which is the typical case if we want to compute tail index
        # Arpack SM is much faster, but it does not always give the right eigenvector (either because LR ≠ SM (happens when the eigenvalue is very positive)
        # Even when it gives the right eigenvalue, it can return a complex eigenvector
        if eigenvector ∈ (:right, :both)
            vals, vecs, info = KrylovKit.eigsolve(𝔸, r0, 1, :LR, maxiter = size(𝔸, 1))
            info.converged == 0 &&  @warn "KrylovKit did not converge"
            η = vals[1]
            f = vecs[1]
        end
        if eigenvector ∈ (:left, :both)
            vals, vecs, info = KrylovKit.eigsolve(adjoint(𝔸), 1, :LR, maxiter = size(𝔸, 1))
            info.converged == 0 &&  @warn "KrylovKit did not converge"
            η = vals[1]
            g = vecs[1]
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
