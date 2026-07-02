"""
    stationary_distribution(𝕋; δ = 0.0, ψ = Ones(size(𝕋, 1)))

Computes the stationary distribution corresponding to the generator matrix `𝕋`.
"""
function stationary_distribution(𝕋::AbstractMatrix; δ = 0.0, ψ = Ones(size(𝕋, 1)))
    size(𝕋, 1) == size(𝕋, 2) || throw(DimensionMismatch("𝕋 must be a square generator matrix"))
    δ >= 0 || throw(ArgumentError("δ needs to be nonnegative"))
    n = size(𝕋, 1)
    ψ = vec(ψ)
    length(ψ) == n || throw(DimensionMismatch("𝕋 and ψ should have the same length"))
    if δ > 0
        g = abs.((δ * I - 𝕋') \ (δ * collect(ψ)))
    else
        η, g = principal_eigenvalue(𝕋')
        abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
    end
    total_mass = sum(g)
    isfinite(total_mass) && total_mass > 0 || throw(ArgumentError("stationary distribution has zero or non-finite mass; pass a positive ψ when δ > 0"))
    g ./ total_mass
end

function _reshape_state_output(X::ContinuousTimeMarkovProcess, x::AbstractVector)
    length(x) == length(X) ||
        throw(DimensionMismatch("state vector has length $(length(x)) but process has length $(length(X))"))
    length(size(X)) == 1 && return x
    return reshape(x, size(X))
end

function _reshape_state_time_output(X::ContinuousTimeMarkovProcess, x::AbstractMatrix)
    size(x, 1) == length(X) ||
        throw(DimensionMismatch("state-time array has $(size(x, 1)) rows but process has length $(length(X))"))
    length(size(X)) == 1 && return x
    return reshape(x, size(X)..., size(x, 2))
end

"""
 computes the stationary distribution corresponding to the continuous-time Markov process X
"""
function stationary_distribution(X::ContinuousTimeMarkovProcess; δ = 0.0, ψ = Ones(length(X)))
    _reshape_state_output(X, stationary_distribution(generator(X); δ = δ, ψ = vec(ψ)))
end
