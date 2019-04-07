module InfinitesimalGenerators
using LinearAlgebra, BandedMatrices, KrylovKit


#========================================================================================

Define Type

========================================================================================#

#struct InfinitesimalGenerator{T, CONTAINER, RAXIS} <: BandedMatrices.AbstractBandedMatrix{T}
#    B::BandedMatrix{T, CONTAINER, RAXIS}
#end
#Lazy.@forward InfinitesimalGenerator.B Base.axes, Base.size, Base.getindex, Base.setindex!, Base.copy
#Base.convert(::Type{T}, 𝔸::InfinitesimalGenerator) where {T <: BandedMatrix}= convert(T, 𝔸.B)
#import Base.+
#(+)(x::InfinitesimalGenerator, y::InfinitesimalGenerator) =  InfinitesimalGenerator(x.B + y.B)
#
#
#Lazy.@forward InfinitesimalGenerator.B LinearAlgebra.svdvals!, LinearAlgebra.factorize
#Lazy.@forward InfinitesimalGenerator.B SparseArrays.sparse
#Lazy.@forward InfinitesimalGenerator.B BandedMatrices.bandeddata, BandedMatrices.bandwidths, BandedMatrices.#data_colrange, BandedMatrices.data_rowrange,  BandedMatrices.MemoryLayout
#@inline BandedMatrices.inbands_getindex(𝔸::InfinitesimalGenerator, u::Integer, k::Integer, j::Integer) = #BandedMatrices.inbands_getindex(𝔸.B, u, k, j)
#@inline BandedMatrices.inbands_getindex(𝔸::InfinitesimalGenerator, k::Integer, j::Integer) = BandedMatrices.inbands_getindex(𝔸.B, k, j)



#========================================================================================

Compute the operator
𝔸f = v_0 * f + v1 * ∂(f) + v2 * ∂∂(f)
Note that
𝔸'g = v_0 * g - ∂(v1 * g) + ∂∂(v2 * g)

========================================================================================#

function operator(x::AbstractVector, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    operator!(𝔸, make_Δ(x), v0, v1, v2)
end

function operator!(𝔸::AbstractMatrix, Δ, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    x, invΔx, invΔxm, invΔxp = Δ
    n = length(x)
    fill!(𝔸, 0.0)
    # construct matrix T. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    for i in 1:n
        if v1[i] >= 0
            𝔸[i, min(i + 1, n)] += v1[i] * invΔxp[i]
            𝔸[i, i] -= v1[i] * invΔxp[i]
        else
            𝔸[i, i] += v1[i] * invΔxm[i]
            𝔸[i, max(i - 1, 1)] -= v1[i] * invΔxm[i]
        end
        𝔸[i, max(i - 1, 1)] += v2[i] * invΔxm[i] * invΔx[i]
        𝔸[i, i] -= v2[i] * 2 * invΔxm[i] * invΔxp[i]
        𝔸[i, min(i + 1, n)] += v2[i] * invΔxp[i] * invΔx[i]
    end
    c = sum(𝔸, dims = 2)
    for i in 1:n
        𝔸[i, i] += v0[i] - c[i]
    end
    return 𝔸
end

function make_Δ(x)
    n = length(x)
    Δxm = zero(x)
    Δxm[1] = x[2] - x[1]
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zero(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δxp[end] = x[n] - x[n-1]
    Δx = (Δxm .+ Δxp) / 2
    return x, 1 ./ Δx, 1 ./ Δxm, 1 ./ Δxp
end

#========================================================================================

Compute the principal eigenvector and eigenvalue of 𝔸

========================================================================================#
function principal_eigenvalue(𝔸::AbstractMatrix; method = :krylov, eigenvector = :right)
    η = nothing
    if method == :krylov
        g, η, f = principal_eigenvalue_krylov(𝔸; eigenvector = eigenvector)
        if η == nothing
            @warn "Krylov Methods Failed"
        end
    end
    if η == nothing
        g, η, f = principal_eigenvalue_BLAS(convert(Matrix{Float64}, 𝔸); eigenvector = eigenvector)
    end
    return clean_eigenvector_left(g), clean_eigenvalue(η), clean_eigenvector_right(f)
end

# I could also use Arpack.eigs but it seems slower
function principal_eigenvalue_krylov(𝔸::AbstractMatrix; eigenvector = :right)
    g, η, f = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        vals, vecs, info = KrylovKit.eigsolve(𝔸, 1, :LR, maxiter = size(𝔸, 1))
        if info.converged > 0
            η = vals[1]
            f = vecs[1]
        end
    end
    if eigenvector ∈ (:left, :both)
        vals, vecs, info = KrylovKit.eigsolve(adjoint(𝔸), 1, :LR, maxiter = size(𝔸, 1))
        if info.converged > 0
            η = vals[1]
            g = vecs[1]
        end
    end 
    return g, η, f
end

function principal_eigenvalue_BLAS(𝔸::AbstractMatrix; eigenvector = :right)
    g, η, f = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        e = eigen(𝔸)
        _, out = findmax(real.(e.values))
        η = e.values[out]
        f = e.vectors[:, out]
    end
    if eigenvector ∈ (:left, :both)
        e = eigen(copy(adjoint(𝔸)))
        _, out = findmax(real.(e.values))
        η = e.values[out]
        g = e.vectors[:, out]
    end 
    return g, η, f
end

clean_eigenvalue(η::Union{Nothing, Real}) = η
function clean_eigenvalue(η::Complex)
    if abs(imag(η) .>= eps())
        @warn "Principal Eigenvalue has some imaginary part $(η)"
    end
    real(η)
end
clean_eigenvector_left(::Nothing) = nothing
clean_eigenvector_left(g::Vector) = abs.(g) ./ sum(abs.(g))
clean_eigenvector_right(::Nothing) = nothing
clean_eigenvector_right(f::Vector) = abs.(f)

#========================================================================================
Solve the PDE backward in time
u(x, T) = ψ(x)
0 = u_t + 𝔸u_t - V(x, t)u +  f(x, t)

using an implicit finite difference scheme, that is
u_T = ψ
u_t = (I - 𝔸dt) \ (u_{t+1} + f dt)
========================================================================================#

function feynman_kac_backward(𝔸::AbstractMatrix; 
	t::AbstractVector = range(0, 100, step = 1/12), 
	ψ::AbstractVector = ones(size(𝔸, 1)), 
	f::Union{AbstractVector, AbstractMatrix} = zeros(size(𝔸, 1)), 
	V::Union{AbstractVector, AbstractMatrix} = zeros(size(𝔸, 1)))
    u = zeros(size(𝔸, 1), length(t))
    u[:, length(t)] = ψ
    if isa(f, AbstractVector) && isa(V, AbstractVector)
        if isa(t, AbstractRange)
            dt = step(t)
            𝔹 = factorize(I + Diagonal(V) * dt - 𝔸 * dt)
            for i in (length(t)-1):(-1):1
                ψ = ldiv!(𝔹, u[:, i+1] .+ f .* dt)
                u[:, i] .= ψ
            end
        else
            for i in (length(t)-1):(-1):1
                dt = t[i+1] - t[i]
                𝔹 = I + Diagonal(V) * dt - 𝔸 * dt
                ψ = 𝔹 \  (u[:, i+1] .+ f .* dt)
                u[:, i] .= ψ
            end
        end
    elseif isa(f, AbstractMatrix) && isa(V, AbstractMatrix)
        for i in (length(t)-1):(-1):1
            dt = t[i+1] - t[i]
            𝔹 = I + Diagonal(V[:, i]) * dt - 𝔸 * dt
            ψ = 𝔹 \ (u[:, i+1] .+ f[:, i] .* dt)
            u[:, i] .= ψ
        end
    else
        error("f and V must be Vectors or Matrices")
    end
    return u
end

#========================================================================================
Solve the PDE forward in time
u(x, 0) = ψ(x)
u_t = 𝔸u - V(x)u + f(x)

using implicit finite difference scheme, that is
u_0 = ψ
u_t = (I - 𝔸dt) \ (u_{t+1} + f dt)
========================================================================================#

function feynman_kac_forward(𝔸::AbstractMatrix; 
	t::AbstractVector = range(0, 100, step = 1/12), 
	ψ::AbstractVector = ones(size(𝔸, 1)), 
	f::AbstractVector = zeros(size(𝔸, 1)), 
	V::AbstractVector = zeros(size(𝔸, 1)))
    u = feynman_kac_backward(𝔸; ψ = ψ, t = - reverse(t), f = f, V = V)
    return u[:,end:-1:1]
end

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
    if abs(η) >= 1e-5
        @warn "Principal Eigenvalue does not seem to be zero"
    end
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
function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector)
    operator(x, μM, σM .* σx .+ μx, 0.5 * σx.^2)
end

# Compute Hansen Scheinkmann decomposition M_t= e^{ηt}f(x_t)W_t
function hansen_scheinkman(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector)
	principal_eigenvalue(generator(x, μx, σx, μM, σM); eigenvector = :right)[2:3]
end

# Compute E[M_t ψ(x_t)|x_0 = x]
function feynman_kac_forward(x::AbstractVector, μx::AbstractVector, σx::AbstractVector,  μM::AbstractVector, σM::AbstractVector; kwargs...)
    feynman_kac_forward(generator(x, μx, σx, μM, σM); kwargs...)
end

##############################################################################
##
## Exported methods and types 
##
##############################################################################
export generator,
principal_eigenvalue,
feynman_kac_backward,
feynman_kac_forward,
stationary_distribution,
hansen_scheinkman
end