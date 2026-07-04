"""
    feynman_kac(𝔸::AbstractMatrix, ts; f, ψ, v, direction = :backward)

Solve the Feynman–Kac PDE associated with a generator matrix `𝔸` on the increasing time
grid `ts`, using implicit Euler time steps.

With `direction = :backward`, returns the solution of

    u(x, ts[end]) = ψ(x)
    0 = uₜ + 𝔸u - v(x, t)u + f(x, t)

or, equivalently, in integral form,

    u(x, t) = E[∫ₜᵀ e^{-∫ₜˢ v(x_u) du} f(x_s) ds + e^{-∫ₜᵀ v(x_u) du} ψ(x_T) | xₜ = x].

With `direction = :forward`, returns the solution of

    u(x, ts[1]) = ψ(x)
    uₜ = 𝔸u - v(x, t)u + f(x, t)

or, equivalently, in integral form,

    u(x, t) = E[∫₀ᵗ e^{-∫₀ˢ v(x_u) du} f(x_s) ds + e^{-∫₀ᵗ v(x_u) du} ψ(xₜ) | x₀ = x].

The inputs `f`, `ψ`, and `v` are vectors of length `size(𝔸, 1)`. The inputs `f`
and `v` may also be matrices with `length(ts)` columns (one per date) or with a
single column (held constant over time).
"""
function feynman_kac(𝔸::AbstractMatrix, ts;
    f::Union{AbstractVector, AbstractMatrix} = zeros(eltype(𝔸), size(𝔸, 1)),
    ψ::AbstractVector = zeros(eltype(𝔸), size(𝔸, 1)),
    v::Union{AbstractVector, AbstractMatrix} = zeros(eltype(𝔸), size(𝔸, 1)),
    direction = :backward)
    size(𝔸, 1) == size(𝔸, 2) || throw(DimensionMismatch("𝔸 must be square matrix"))
    size(𝔸, 1) == size(f, 1) || throw(DimensionMismatch("𝔸 and f should have the same number of rows"))
    size(𝔸, 1) == length(ψ) || throw(DimensionMismatch("𝔸 and ψ should have the same number of rows"))
    size(𝔸, 1) == size(v, 1) || throw(DimensionMismatch("𝔸 and v should have the same number of rows"))
    size(f, 2) ∈ (1, length(ts)) ||  throw(DimensionMismatch("The number of columns in f should equal the length of ts"))
    size(v, 2) ∈ (1, length(ts)) ||  throw(DimensionMismatch("The number of columns in v should equal the length of ts"))
    direction ∈ (:forward, :backward) || throw(ArgumentError("Direction must be :backward or :forward"))
    issorted(ts) || throw(ArgumentError("`ts` must be increasing"))
    if direction == :forward
        f_reverse = ndims(f) == 2 ? @view(f[:, end:-1:1]) : f
        v_reverse = ndims(v) == 2 ? @view(v[:, end:-1:1]) : v
        u = feynman_kac(𝔸, - reverse(ts); ψ = ψ, f = f_reverse, v = v_reverse, direction = :backward)
        return u[:,end:-1:1]
    end
    # direction is backward
    # one-column f/v are held constant over time; otherwise column i is the value at ts[i]
    f_col = i -> ndims(f) == 1 ? f : view(f, :, size(f, 2) == 1 ? 1 : i)
    v_col = i -> ndims(v) == 1 ? v : view(v, :, size(v, 2) == 1 ? 1 : i)
    T = float(promote_type(eltype(𝔸), eltype(f), eltype(ψ), eltype(v), eltype(ts)))
    u = zeros(T, size(𝔸, 1), length(ts))
    u[:, end] = ψ
    if isa(ts, AbstractRange) && (ndims(v) == 1 || size(v, 2) == 1)
        # constant time step and time-invariant v: factorize once
        dt = step(ts)
        B = factorize(I + (Diagonal(v_col(1)) - 𝔸) * dt)
        rhs = Vector{T}(undef, size(𝔸, 1))
        # Use in-place solves for tridiagonal/banded factors, but keep
        # a fallback for sparse factors such as CHOLMOD without ldiv!.
        can_ldiv = hasmethod(ldiv!, Tuple{typeof(B), typeof(rhs)})
        for i in (length(ts)-1):(-1):1
            @views copyto!(rhs, u[:, i + 1])
            rhs .+= f_col(i) .* dt
            if can_ldiv
                ldiv!(B, rhs)
                @views copyto!(u[:, i], rhs)
            else
                @views u[:, i] .= B \ rhs
            end
        end
    else
        # non-constant time step or time-varying v: refactorize at each step
        for i in (length(ts)-1):(-1):1
            dt = ts[i+1] - ts[i]
            B = I + (Diagonal(v_col(i)) - 𝔸) * dt
            u[:, i] = B \ (u[:, i+1] .+ f_col(i) .* dt)
        end
    end
    return u
end
