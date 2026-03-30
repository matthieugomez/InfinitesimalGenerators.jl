"""
Compute the principal eigenvalue and eigenvector of a Metzler matrix 𝕋
(i.e. a matrix with non-negative off-diagonal entries).

By Perron-Frobenius, the eigenvalue η with largest real part is real,
and the corresponding eigenvector r is strictly positive.

Two cases:
1. If rows or columns sum to zero (𝕋 is a generator), then η = 0.
   The eigenvector is found by solving 𝕋r = 0 with r[1] = 1.
2. Otherwise, η is found by inverse iteration with Rayleigh quotient updates.
   At each step, we solve (𝕋 - σI)w = r, normalize w, and update the shift σ
   via the Rayleigh quotient σ = r'𝕋r. This converges cubically to η.
   The initial shift is the max row sum — a Gershgorin upper bound on η —
   which guarantees that η is the eigenvalue closest to the shift,
   so inverse iteration converges to the right eigenvalue.
   For tridiagonal 𝕋, each iteration costs O(n).
"""
function principal_eigenvalue(𝕋; r0 = ones(size(𝕋, 1)), η0 = nothing, maxiter = 100, tol = 1e-12)
    if (maximum(abs.(sum(𝕋, dims = 1))) < 1e-9) || (maximum(abs.(sum(𝕋, dims = 2))) < 1e-9)
        # rows or columns sum to zero → η = 0, solve 𝕋r = 0 with r[1] = 1
        if 𝕋 isa Tridiagonal
            r = [1.0 ; - Tridiagonal(𝕋.dl[2:end], 𝕋.d[2:end], 𝕋.du[2:end]) \ vec(𝕋[2:end, 1])]
        else
            r = [1.0 ; - 𝕋[2:end, 2:end] \ vec(𝕋[2:end, 1])]
        end
        return 0.0, abs.(r)
    else
        # Inverse iteration with Rayleigh quotient updates
        r = collect(float.(r0))
        r ./= sqrt(r' * r)
        if η0 !== nothing
            η = float(η0)
        else
            η = maximum(sum(𝕋, dims = 2))
        end
        for _ in 1:maxiter
            w = (𝕋 - η * I) \ r
            r = w ./ sqrt(w' * w)
            η_new = r' * (𝕋 * r)
            if abs(η_new - η) < tol * (1 + abs(η_new))
                return η_new, abs.(r)
            end
            η = η_new
        end
        @warn "Inverse iteration did not converge in $maxiter iterations"
        return η, abs.(r)
    end
end
