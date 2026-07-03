function _resolve_rebirth_argument(rebirth, ψ, default)
    if ψ !== nothing
        Base.depwarn("the `ψ` keyword argument of `stationary_distribution` is deprecated; use `rebirth` instead", :stationary_distribution)
        rebirth === nothing || throw(ArgumentError("pass only one of `rebirth` and `ψ`"))
        return ψ
    end
    return rebirth === nothing ? default : rebirth
end

"""
    stationary_distribution(𝔸::AbstractMatrix; δ = 0.0, rebirth = Ones(size(𝔸, 1)))

Compute the stationary distribution corresponding to the generator matrix `𝔸` — the
solution of the Kolmogorov forward equation `𝔸'g = 0`, normalized to sum to one.

With a death rate `δ > 0` and a rebirth distribution `ψ` given by `rebirth`, returns
instead the resolvent `(δI - 𝔸')⁻¹ δψ` — the stationary distribution when units die at
rate `δ` and are reborn with distribution `ψ`.
"""
function stationary_distribution(𝔸::AbstractMatrix; δ = 0.0, rebirth = nothing, ψ = nothing)
    rebirth = _resolve_rebirth_argument(rebirth, ψ, Ones(size(𝔸, 1)))
    size(𝔸, 1) == size(𝔸, 2) || throw(DimensionMismatch("𝔸 must be a square generator matrix"))
    δ >= 0 || throw(ArgumentError("δ needs to be nonnegative"))
    n = size(𝔸, 1)
    rebirth = vec(rebirth)
    length(rebirth) == n || throw(DimensionMismatch("𝔸 and rebirth should have the same length"))
    if δ > 0
        g = abs.((δ * I - 𝔸') \ (δ * collect(rebirth)))
    else
        η, g = principal_eigenvalue(𝔸')
        abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
    end
    total_mass = sum(g)
    isfinite(total_mass) && total_mass > 0 || throw(ArgumentError("stationary distribution has zero or non-finite mass; pass a positive rebirth distribution when δ > 0"))
    g ./ total_mass
end
