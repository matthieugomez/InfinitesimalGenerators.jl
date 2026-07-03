function _default_atol(𝔸::AbstractMatrix)
    T = float(real(eltype(𝔸)))
    sqrt(eps(T)) * max(one(T), maximum(abs, diag(𝔸)))
end

"""
    check_metzler(𝔸::AbstractMatrix; atol = nothing)

Check that `𝔸` is a square Metzler matrix — nonnegative off-diagonal entries, row sums
unrestricted. This is the Perron–Frobenius hypothesis of [`principal_eigenvalue`](@ref):
generators and tilted generators both satisfy it.

A violation beyond `atol` emits a warning reporting its size rather than throwing — a
matrix that is close to Metzler often still yields accurate results. Non-square or empty
matrices throw an `ArgumentError`. Returns `𝔸`. Not exported.
"""
function check_metzler(𝔸::AbstractMatrix; atol = nothing)
    size(𝔸, 1) == size(𝔸, 2) ||
        throw(ArgumentError("matrix must be square"))
    size(𝔸, 1) > 0 ||
        throw(ArgumentError("matrix cannot be empty"))
    if atol === nothing
        atol = _default_atol(𝔸)
    end
    m = _offdiagonal_minimum(𝔸)
    if m < -atol
        @warn "matrix is not Metzler: minimum off-diagonal entry is $m (tolerance $atol). " *
            "Perron–Frobenius does not strictly apply, but results may still be accurate if the violation is small."
    end
    return 𝔸
end

"""
    check_generator(𝔸::AbstractMatrix; atol = nothing)

Check that `𝔸` is a valid generator (transition-rate) matrix: square, with nonnegative
off-diagonal entries and rows summing to zero — that is, a Metzler matrix whose rows sum
to zero.

A violation beyond `atol` emits a warning reporting its size rather than throwing — a
matrix that is close to a generator often still yields accurate results. Non-square or
empty matrices throw an `ArgumentError`. Returns `𝔸`, so the call can be used inline.

`atol` defaults to `sqrt(eps) * max(1, maximum(abs, diag(𝔸)))` (for a valid generator, the
largest entry in absolute value is on the diagonal). The off-diagonal check is
`minimum(𝔸 - Diagonal(diag(𝔸)))`, which preserves the structure of the matrix — sparse
matrices are checked in one pass over their stored entries, tridiagonal ones over their
bands.

Note that a *tilted* generator ([`tilted_generator`](@ref)`(m, ξ)`) is intentionally not a
generator in this sense: it has nonnegative off-diagonal entries but nonzero row sums —
that is why `E[e^{ξ mₜ}]` grows at rate `Λ(ξ)` instead of remaining a probability.
"""
function check_generator(𝔸::AbstractMatrix; atol = nothing)
    if atol === nothing
        atol = _default_atol(𝔸)
    end
    check_metzler(𝔸; atol = atol)
    T = float(real(eltype(𝔸)))
    r = maximum(abs, 𝔸 * ones(T, size(𝔸, 2)))
    if r > atol
        @warn "matrix is not a generator (transition-rate) matrix: maximum |row sum| is $r (tolerance $atol). " *
            "Probability is not conserved, but results may still be accurate if the violation is small."
    end
    return 𝔸
end

_offdiagonal_minimum(𝔸::AbstractMatrix) = minimum(𝔸 - Diagonal(diag(𝔸)))

# generic reductions over banded types visit all n² entries; use the bands directly
_offdiagonal_minimum(𝔸::Tridiagonal) =
    isempty(𝔸.dl) ? zero(eltype(𝔸)) : min(minimum(𝔸.dl), minimum(𝔸.du))

_offdiagonal_minimum(𝔸::Diagonal) = zero(eltype(𝔸))

# the off-diagonal minimum of a real matrix is transpose-invariant
_offdiagonal_minimum(𝔸::Union{Adjoint{<:Real}, Transpose{<:Real}}) = _offdiagonal_minimum(parent(𝔸))
