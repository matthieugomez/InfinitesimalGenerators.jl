module InfinitesimalGenerators
using LinearAlgebra, Arpack, Roots

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

# Compute generator 𝔸f = E[d(Mf(x))]
function generator(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number})
    operator(x, μM, μx .+ σM .* σx, 0.5 * σx.^2)
end

# Compute Hansen Scheinkmann decomposition M_t= e^{ηt}f(x_t)\hat{M}_t
function hansen_scheinkman(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; eigenvector = :right, symmetrize = false)
    principal_eigenvalue(generator(x, μx, σx, μM, σM); which = :LR, eigenvector = eigenvector)
end

# Compute E[M_t ψ(x_t)|x_0 = x]
function feynman_kac_forward(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number},  μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; kwargs...)
    feynman_kac_forward(generator(x, μx, σx, μM, σM); kwargs...)
end


##############################################################################
##
## Tail Index
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
    ζ = find_zero(mgf_longrun(x, μx, σx, μM, σM; δ = δ, ρ = ρ), (1e-3, 40.0))
    return ζ
end

# ξ -> 𝔸(ξ)
function mgf_generator(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; δ::Number = 0.0,  ρ::Number = 0.0)
    ξ -> operator(x, ξ .* μM .+ 0.5 * ξ * (ξ - 1) .* σM.^2 .- δ,  μx .+ ξ .* σM .* ρ .* σx, 0.5 * σx.^2)
end

# ξ -> lim(log(E[M_t^ξ|x_0 = x])/t)
function mgf_longrun(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; δ::Number = 0.0,  ρ::Number = 0.0)
    ξ -> begin
        M = mgf_generator(x, μx, σx, μM, σM; δ = δ, ρ = ρ)(ξ)
        out = nothing
        while out === nothing
            try 
                out = principal_eigenvalue(M; which = :LR, eigenvector = :right)[2]
            catch
                # LR fails when the LR eigenvalue is very close to zero, i.e. around the solution I'm interested in it. 
                # in this case, the SM eigenvalue is the LR eigenvalue, so I can just use SM
                # however, i cannot do SM everywhere: for large value of ξ, we have LR > 0 and there may be one close to zero.
                # This fix does not work if LR fails in a region where LR ≠ SM
                # in this case i should just restart LR until it works
                out = principal_eigenvalue(M; which = :SM, eigenvector = :right)[2]
            end
        end
        return out
    end
end

# Compute first derivative of ξ -> lim(log(E[M_t^ξ|x_0 = x])/t)
function ∂mgf_longrun(x::AbstractVector{<:Number}, μx::AbstractVector{<:Number}, σx::AbstractVector{<:Number}, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; δ::Number = 0.0,  ρ::Number = 0.0)
    return ξ -> begin
        g, η, f = principal_eigenvalue(mgf_generator(x, μx, σx, μM, σM; δ = δ, ρ = ρ)(ξ); which = :LR, eigenvector = :both)
        ∂𝔸 = operator(x, μM .+ (η - 1/2) .* σM.^2, σM .* ρ .* σx, zeros(length(x)))
        (g' * ∂𝔸 * f) / (g' * f)
    end
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
mgf_generator,
mgf_longrun,
hansen_scheinkman,
tail_index
end