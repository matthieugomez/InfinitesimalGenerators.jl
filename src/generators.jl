
##############################################################################
##
## Markov Process
##
## must define generator!(x)
## which corresponds to 𝔸f = E[df]
##
##############################################################################

abstract type MarkovProcess end
# must define generator!

function generator(x::MarkovProcess)
    deepcopy(generator!(x))
end

function stationary_distribution(x::MarkovProcess)
    g, η, _ = principal_eigenvalue(generator!(x); which = :SM, eigenvector = :left)
    abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
    return g
end

# Death rate δ and reinjection ψ
function stationary_distribution(x::MarkovProcess, δ::Number, ψ::AbstractVector{<:Number})
    clean_eigenvector_left((δ * I - generator!(x)') \ (δ * ψ))
end

"""
If direction = :backward
compute `u(x, t) = E[∫t^T e^{-∫ts V(x_τ, τ)dτ}f(x_s, s)ds + e^{-∫tT V(x_τ, τ)dτ}ψ(x_T)|x_t = x]`
If direction = :forward
compute `u(x, t)= E[∫0^t e^{-∫0^s V(x_τ)dτ}f(x_s)ds + e^{-∫0^tV(x_τ)dτ}ψ(x_t)|x_0 = x]`
"""
function feynman_kac(x::MarkovProcess; kwargs...)
    feynman_kac(generator!(x); kwargs...)
end

##############################################################################
##
## Multiplicative Functional
##
## must define generator!(M, ξ::Real)
## which corresponds to 𝔸f = E[d(M^ξf)]
## Must also define length
##############################################################################


abstract type MultiplicativeFunctional end
 
function generator(M::MultiplicativeFunctional, ξ = 1.0)
    deepcopy(generator!(M, ξ))
end

# ξ -> lim log(E[M^\xi]) / t
function cgf_longrun(M::MultiplicativeFunctional; which = :LR, eigenvector = :right, r0 = Ones(length(M)))
    ξ -> principal_eigenvalue(generator!(M, ξ); which = which, eigenvector = eigenvector, r0 = r0)
end

# Compute Hansen Scheinkmann decomposition M_t= e^{ηt}f(x_t)\hat{M}_t
function hansen_scheinkman(M::MultiplicativeFunctional; which = :LR, eigenvector = :right)
    cgf_longrun(M, eigenvector = eigenvector)(1.0)
end

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
function tail_index(M::MultiplicativeFunctional; kwargs...)
    find_root(ξ -> generator!(M, ξ); kwargs...)
end

""" 
If direction = :forward
compute `E[M_t ψ(x_t)|x_0 = x]`
"""
function feynman_kac(M::MultiplicativeFunctional; kwargs...)
    feynman_kac(generator!(M); kwargs...)
end
