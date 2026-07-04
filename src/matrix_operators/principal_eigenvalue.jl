# exceptions thrown by the various `\` methods (dense LAPACK, tridiagonal, sparse UMFPACK)
# when the matrix is exactly singular; ZeroPivotException only exists on Julia ≥ 1.7
_is_singular_exception(e) =
    e isa SingularException || e isa LAPACKException ||
    (isdefined(LinearAlgebra, :ZeroPivotException) && e isa LinearAlgebra.ZeroPivotException)

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
1. If rows or columns sum to zero (𝔸 is a generator, or its transpose), then η = 0.
   The eigenvector is found by solving 𝔸r = 0 with r[1] = 1; the row-sum test is scaled
   by the size of the diagonal, so it is robust to the rounding of generator assembly.
   If that solve is singular, the generator is reducible (it has states that do not
   communicate), the stationary eigenvector is not unique, and an `ArgumentError` says so.
2. Otherwise, η is found by inverse iteration with Rayleigh quotient updates.
   At each step, we solve (𝔸 - σI)w = r, normalize w, and update the shift σ
   via the Rayleigh quotient σ = r'𝔸r. This converges cubically to a nearby
   eigenvalue. The initial shift is the max row sum — a Gershgorin upper bound
   on the real part of the spectrum. When 𝔸 has a real spectrum (e.g. the
   tridiagonal generators of 1-D diffusions, which are similar to symmetric
   matrices), the eigenvalue nearest that shift is the principal one, so the
   iteration lands on η. For a general Metzler matrix with complex eigenvalues
   this is not guaranteed; pass `η0` to start from a known bound if needed.
   If a shift lands exactly on an eigenvalue (e.g. the Gershgorin bound when all
   row sums are equal, so the bound is attained), it is nudged by `tol` to keep
   the solve nonsingular; the Rayleigh quotient then recovers the exact eigenvalue.
   For tridiagonal 𝔸, each iteration costs O(n).
"""
function principal_eigenvalue(𝔸::AbstractMatrix; r0 = ones(size(𝔸, 1)), η0 = nothing, maxiter = 100, tol = 1e-12)
    check_metzler(𝔸)
    o = float(real(one(eltype(𝔸))))
    # scaled like `_default_atol` but tighter: it only needs to absorb the rounding of
    # generator assembly, without misclassifying weakly tilted generators as generators
    atol_generator = 1000 * eps(o) * max(o, maximum(abs, diag(𝔸)))
    if (maximum(abs.(sum(𝔸, dims = 1))) < atol_generator) || (maximum(abs.(sum(𝔸, dims = 2))) < atol_generator)
        # rows or columns sum to zero (a generator or its transpose) → η = 0,
        # solve 𝔸r = 0 with r[1] = 1
        r = try
            if 𝔸 isa Tridiagonal
                [1.0 ; - Tridiagonal(𝔸.dl[2:end], 𝔸.d[2:end], 𝔸.du[2:end]) \ vec(𝔸[2:end, 1])]
            else
                [1.0 ; - 𝔸[2:end, 2:end] \ collect(𝔸[2:end, 1])]
            end
        catch e
            _is_singular_exception(e) || rethrow()
            throw(ArgumentError("could not compute the zero eigenvector: the generator appears " *
                "reducible (it has states that do not communicate), so the stationary eigenvector " *
                "is not unique. Restrict the process to a single communicating class."))
        end
        return zero(o), abs.(r)
    else
        # Inverse iteration with Rayleigh quotient updates
        r = collect(float.(r0))
        r ./= sqrt(r' * r)
        η = η0 !== nothing ? float(η0) : maximum(sum(𝔸, dims = 2))
        for _ in 1:maxiter
            w = nothing
            try
                w0 = (𝔸 - η * I) \ r
                all(isfinite, w0) && (w = w0)
            catch e
                _is_singular_exception(e) || rethrow()
            end
            if w === nothing
                # η is exactly an eigenvalue: nudge the shift so the solve is nonsingular;
                # the Rayleigh quotient below recovers the exact eigenvalue
                w = (𝔸 - (η + tol * (1 + abs(η))) * I) \ r
            end
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
