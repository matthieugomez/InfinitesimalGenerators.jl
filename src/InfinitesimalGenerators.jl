module InfinitesimalGenerators
using LinearAlgebra, Arpack, KrylovKit, Roots

include("operators.jl")
include("feynman_kac.jl")

#========================================================================================

For a Markov Process x:
dx = μx dt + σx dZ_t

========================================================================================#

# Compute generator 𝔸f = E[df(x)]
function generator(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number})
    operator(x, zeros(length(x)), μx, 0.5 * σx.^2)
end

# Stationary Distribution of x
function stationary_distribution(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number})
    g, η, _ = principal_eigenvalue(generator(x, μx, σx); which = :SM, eigenvector = :left)
    abs(η) >= 1e-5 && @warn "Principal Eigenvalue does not seem to be zero"
    return g
end

# Stationary Distribution of x with death rate δ and reinjection ψ
function stationary_distribution(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, δ::Number, ψ::AbstractVector{<:Number})
    clean_eigenvector_left((δ * I - adjoint(generator(x, μx, σx))) \ (δ * ψ))
end

# Compute u(x_t, t) = E[∫t^T e^{-∫ts V(x_τ, τ)dτ}f(x_s, s)ds + e^{-∫tT V(x_τ, τ)dτ}ψ(x_T)|x_t = x]
function feynman_kac_backward(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}; kwargs...)
    feynman_kac_backward(generator(x, μx, σx); kwargs...)
end

# Compute u(x, t)= E[∫0^t e^{-∫0^s V(x_τ)dτ}f(x_s)ds + e^{-∫0^tV(x_τ)dτ} ψ(x_t)|x_0 = x]
function feynman_kac_forward(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}; kwargs...)
    feynman_kac_forward(generator(x, μx, σx); kwargs...)
end

#========================================================================================

For a Markov Process x:
dx = μx dt + σx dZt
and a multiplicative functional M:
dM/M = μM dt + σM dZt

========================================================================================#



##############################################################################
##
## MultiiplicativeFunctionals
##
##############################################################################
struct MultiplicativeFunctional
    𝔸::Tridiagonal
    Δ::Tuple{<:AbstractVector, <:AbstractVector, <:AbstractVector, <:AbstractVector}
    μx::AbstractVector{<:Number}
    σx::AbstractVector{<:Number}
    μM::AbstractVector{<:Number}
    σM::AbstractVector{<:Number}
    δ::Number
    ρ::Number
end

function MultiplicativeFunctional(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; δ::Number = 0.0,  ρ::Number = 0.0)
    n = length(x)
    𝔸 = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    Δ = make_Δ(x)
    MultiplicativeFunctional(𝔸, Δ, μx, σx, μM, σM, δ, ρ)
end

# ξ -> 𝔸(ξ)
function generator(M::MultiplicativeFunctional)
    operator!(M.𝔸, M.Δ, M.μM .- M.δ,  M.μx .+ M.σM .* M.ρ .* M.σx, 0.5 * M.σx.^2)
end

# Compute Hansen Scheinkmann decomposition M_t= e^{ηt}f(x_t)\hat{M}_t
function hansen_scheinkman(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; δ::Number = 0.0,  ρ::Number = 0.0, eigenvector = :right)
    hansen_scheinkman(MultiplicativeFunctional(x, μx, σx, μM, σM; δ = δ, ρ = ρ), eigenvector = eigenvector)
end
function hansen_scheinkman(M::MultiplicativeFunctional; eigenvector = :right)
    principal_eigenvalue(generator(M); which = :LR, eigenvector = eigenvector)
end

# Compute E[M_t ψ(x_t)|x_0 = x]
function feynman_kac_forward(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; δ::Number = 0.0,  ρ::Number = 0.0, kwargs...)
    feynman_kac_forward(MultiplicativeFunctional(x, μx, σx, μM, σM; δ = δ, ρ = ρ); kwargs...)
end

function feynman_kac_forward(M::MultiplicativeFunctional; kwargs...)
    feynman_kac_forward(generator(M); kwargs...)
end

##############################################################################
##
## CGF
##
##############################################################################

function generator(M::MultiplicativeFunctional, ξ)
    operator!(M.𝔸, M.Δ, ξ .* M.μM .+ 0.5 * ξ * (ξ - 1) .* M.σM.^2 .- M.δ,  M.μx .+ ξ .* M.σM .* M.ρ .* M.σx, 0.5 * M.σx.^2)
end

# ξ -> \lim log(E[M^\xi]) / t
function cgf_longrun(M::MultiplicativeFunctional, ξ; eigenvector = :right)
    principal_eigenvalue(generator(M, ξ), which = :LR, eigenvector = eigenvector)
end

# Compute first derivative of ξ -> lim(log(E[M_t^ξ|x_0 = x])/t)
function ∂cgf_longrun(M::MultiplicativeFunctional, ξ::Number)
    g, η, f = principal_eigenvalue(generator(M, ξ); which = :LR, eigenvector = :both)
    ∂𝔸 = operator(x, μM .+ (η - 1/2) .* σM.^2, σM .* ρ .* σx, zeros(length(x)))
    return (g' * ∂𝔸 * f) / (g' * f)
end

##############################################################################
##
## Tail Indices
##
##############################################################################

# Compute tail index of the process M given by
# dM/M = μ dt + σ dW_t
# with death rate δ
function tail_index(μ::Number, σ::Number; δ::Number = 0)
    if σ > 0
        (1 - 2 * μ / σ^2 + sqrt((1- 2 * μ / σ^2)^2 + 8 * δ / σ^2)) / 2
    else
        δ / μ
    end
end

# Compute tail index of the process M given by
# dM/M = μM(x) dt + σM(x) dZt
# dx = μx dt + σx dZt
# with death rate δ
function tail_index(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; δ::Number = 0.0,  ρ::Number = 0.0)
    M = MultiplicativeFunctional(x, μx, σx, μM, σM; δ = δ, ρ = ρ)
    find_zero(ξ -> cgf_longrun(M, ξ)[2], (1e-3, 40.0))
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
MultiplicativeFunctional,
hansen_scheinkman,
cgf_longrun,
∂cgf_longrun,
tail_index
end