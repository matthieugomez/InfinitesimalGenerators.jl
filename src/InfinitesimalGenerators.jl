module InfinitesimalGenerators
using LinearAlgebra, BandedMatrices, KrylovKit

include("operators.jl")
include("feynman_kac.jl")

#========================================================================================

For a Markov Process x:
dx = μx dt + σx dZ_t

========================================================================================#

# Compute generator 𝔸f = E[df(x)]
function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)
    operator(x, zeros(length(x)), μx, 0.5 * σx.^2)
end

# Stationary Distribution of x
function stationary_distribution(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)
    g, η, _ = principal_eigenvalue(generator(x, μx, σx); eigenvector = :left)
    abs(η) >= 1e-5 && @warn "Principal Eigenvalue does not seem to be zero"
    return g
end

# Stationary Distribution of x with death rate δ and reinjection ψ
function stationary_distribution(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, δ, ψ)
    clean_eigenvector_left((δ * I - adjoint(generator(x, μx, σx))) \ (δ * ψ))
end

# Compute u(x_t, t) = E[∫t^T e^{-∫ts V(x_τ, τ)dτ}f(x_s, s)ds + e^{-∫tT V(x_τ, τ)dτ}ψ(x_T)|x_t = x]
function feynman_kac_backward(x::AbstractVector, μx::AbstractVector, σx::AbstractVector; kwargs...)
    feynman_kac_backward(generator(x, μx, σx); kwargs...)
end

# Compute u(x, t)= E[∫0^t e^{-∫0^s V(x_τ)dτ}f(x_s)ds + e^{-∫0^tV(x_τ)dτ} ψ(x_t)|x_0 = x]
function feynman_kac_forward(x::AbstractVector, μx::AbstractVector, σx::AbstractVector; kwargs...)
    feynman_kac_forward(generator(x, μx, σx); kwargs...)
end

#========================================================================================

For a Markov Process x:
dx = μx dt + σx dZt
and a multiplicative functional M:
dM/M = μM dt + σM dZt

========================================================================================#

# Compute generator 𝔸f = E[d(Mf(x))]
function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector; symmetrize = false)
    𝔸 = operator(x, μM, σM .* σx .+ μx, 0.5 * σx.^2)
    if symmetrize
        g = stationary_distribution(x, μx, σx)
        𝔸 = Symmetric(Diagonal(sqrt.(g))' * 𝔸 * Diagonal(1 ./ sqrt.(g)))
    end
    return 𝔸
end

# Compute Hansen Scheinkmann decomposition M_t= e^{ηt}f(x_t)W_t
function hansen_scheinkman(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector; symmetrize = false)
    principal_eigenvalue(generator(x, μx, σx, μM, σM; symmetrize = symmetrize); eigenvector = :right)[2:3]
end

# Compute E[M_t ψ(x_t)|x_0 = x]
function feynman_kac_forward(x::AbstractVector, μx::AbstractVector, σx::AbstractVector,  μM::AbstractVector, σM::AbstractVector; kwargs...)
    feynman_kac_forward(generator(x, μx, σx, μM, σM); kwargs...)
end

# Compute tail index of the process M given by
# dM/M = μ dt + σ dW_t
# with death rate δ
function tail_index(μ, σ, δ = 0)
    if σ > 0
        (1 - 2 * μ / σ^2 + sqrt((1- 2 * μ / σ^2)^2 + 8 * δ / σ^2)) / 2
    else
        δ / μ
    end
end



# Compute tail index of the process M given by
# dM/M = μM(x) dt + νM(x) dW_t
# dx = μ(x) dt + σ(x) dW_t
# with death rate δ
function tail_index(x, μx, σx, μM, σM, δ = 0.0)
    ζ = find_zero(ξ -> hansen_scheinkman(x, μx, σx, ξ .* μM .+ 0.5 * ξ * (ξ - 1) .* σM.^2 .- δ, ξ .* σM)[1], (1e-6, 10.0))
    out = hansen_scheinkman(x, μx, σx, ζ .* μM .+ 0.5 * ζ * (ζ - 1) .* σM.^2, ζ .* σM)
    (abs(out) > 1e-3) && @warn "could not find zero power law"
    return ζ
end


##############################################################################
##
## Exported methods and types 
##
##############################################################################

export 
generator,
principal_eigenvalue,
feynman_kac_backward,
feynman_kac_forward,
stationary_distribution,
hansen_scheinkman,
tail_index
end