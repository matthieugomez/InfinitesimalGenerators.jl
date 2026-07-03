"""
    principal_eigenvalue(𝔸; r0 = ones(size(𝔸, 1)), η0 = nothing, maxiter = 100, tol = 1e-12)

Compute the principal eigenvalue and eigenvector of a Metzler matrix 𝔸
(i.e. a matrix with non-negative off-diagonal entries), returned as a tuple `(η, r)`.

By Perron-Frobenius, the eigenvalue η with largest real part is real,
and the corresponding eigenvector r is strictly positive. Metzlerity is what makes this
hold — generators and tilted generators both satisfy it — so the matrix is checked on
entry, with a warning if an off-diagonal entry is negative (results may still be accurate
when the violation is small).

Two cases:
1. If rows or columns sum to zero (𝔸 is a generator), then η = 0.
   The eigenvector is found by solving 𝔸r = 0 with r[1] = 1.
2. Otherwise, η is found by inverse iteration with Rayleigh quotient updates.
   At each step, we solve (𝔸 - σI)w = r, normalize w, and update the shift σ
   via the Rayleigh quotient σ = r'𝔸r. This converges cubically to a nearby
   eigenvalue. The initial shift is the max row sum — a Gershgorin upper bound
   on the real part of the spectrum. When 𝔸 has a real spectrum (e.g. the
   tridiagonal generators of 1-D diffusions, which are similar to symmetric
   matrices), the eigenvalue nearest that shift is the principal one, so the
   iteration lands on η. For a general Metzler matrix with complex eigenvalues
   this is not guaranteed; pass `η0` to start from a known bound if needed.
   For tridiagonal 𝔸, each iteration costs O(n).
"""
function principal_eigenvalue(𝔸::AbstractMatrix; r0 = ones(size(𝔸, 1)), η0 = nothing, maxiter = 100, tol = 1e-12)
    check_metzler(𝔸)
    if (maximum(abs.(sum(𝔸, dims = 1))) < 1e-9) || (maximum(abs.(sum(𝔸, dims = 2))) < 1e-9)
        # rows or columns sum to zero → η = 0, solve 𝔸r = 0 with r[1] = 1
        if 𝔸 isa Tridiagonal
            r = [1.0 ; - Tridiagonal(𝔸.dl[2:end], 𝔸.d[2:end], 𝔸.du[2:end]) \ vec(𝔸[2:end, 1])]
        else
            r = [1.0 ; - 𝔸[2:end, 2:end] \ collect(𝔸[2:end, 1])]
        end
        return 0.0, abs.(r)
    else
        # Inverse iteration with Rayleigh quotient updates
        r = collect(float.(r0))
        r ./= sqrt(r' * r)
        if η0 !== nothing
            η = float(η0)
        else
            η = maximum(sum(𝔸, dims = 2))
        end
        for _ in 1:maxiter
            w = (𝔸 - η * I) \ r
            r = w ./ sqrt(w' * w)
            η_new = r' * (𝔸 * r)
            if abs(η_new - η) < tol * (1 + abs(η_new))
                return η_new, abs.(r)
            end
            η = η_new
        end
        @warn "Inverse iteration did not converge in $maxiter iterations"
        return η, abs.(r)
    end
end
