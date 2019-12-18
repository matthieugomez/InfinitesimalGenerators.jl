
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
function tail_index(M::MultiplicativeFunctional; which = :SM, xatol = 1e-2, verbose = false, r0 = ones(length(M)), kwargs...)
    out = 0.0
    if which == :SM
        try
            # SM is so much faster. So try if it works.
            f = ξ -> begin
                out = cgf_longrun(M; which = :SM, r0 = r0)(ξ)
                eltype(out[3]) <: Float64 && copyto!(r0, out[3])
                verbose && @show (:SM, ξ, out[2])
                return out[2]
            end
            D = ξ -> DiffEqDiffTools.finite_difference_derivative(f, ξ)
            out = find_zero((f, D), 1.0, Roots.Newton(); xatol = xatol, kwargs...)
            out2 = cgf_longrun(M; which = :LR, r0 = r0)(out)[2]
            if abs(out2) > 1e-2 
                @warn "Algorithm looking for SM eigenvalue = 0 converged to ζ = $out. However, the :LR eigenvalue for this ζ is  $out2"
                throw("there is an error")
            end
        catch
            which = :LR
        end
    end
    if which == :LR
        f = ξ -> begin
            out = cgf_longrun(M; which = :LR, r0 = r0)(ξ)
            eltype(out[3]) <: Float64 && copyto!(r0, out[3])
            verbose && @show (:LR, ξ, out[2])
            return out[2]
        end
        D = ξ -> DiffEqDiffTools.finite_difference_derivative(f, ξ)
        try
            out = find_zero((f, D), 1.0, Roots.Newton(); xatol = xatol, kwargs...)
        catch
            out = find_zero((f, D), (1e-2, 10.0); xatol = xatol, kwargs...)
        end
    end
    return out
end

""" 
If direction = :forward
compute `E[M_t ψ(x_t)|x_0 = x]`
"""
function feynman_kac(M::MultiplicativeFunctional; kwargs...)
    feynman_kac(generator!(M); kwargs...)
end
