module InfinitesimalGenerators
using LinearAlgebra, SparseArrays, Lazy, BandedMatrices, KrylovKit

struct InfinitesimalGenerator{T, CONTAINER, RAXIS} <: BandedMatrices.AbstractBandedMatrix{T}
    B::BandedMatrix{T, CONTAINER, RAXIS}
end
Lazy.@forward InfinitesimalGenerator.B Base.axes, Base.size, Base.getindex, Base.setindex!, LinearAlgebra.svdvals!, LinearAlgebra.factorize, SparseArrays.sparse, BandedMatrices.bandeddata, BandedMatrices.bandwidths, BandedMatrices.data_colrange, BandedMatrices.data_rowrange,  BandedMatrices.MemoryLayout, Base.copy

@inline inbands_getindex(𝔸::InfinitesimalGenerator, u::Integer, k::Integer, j::Integer) = BandedMatrices.inbands_getindex(𝔸.B, u, k, j)
@inline inbands_getindex(𝔸::InfinitesimalGenerator, k::Integer, j::Integer) = BandedMatrices.inbands_getindex(𝔸.B, k, j)
Base.convert(::Type{T}, 𝔸::InfinitesimalGenerator) where {T <: BandedMatrix}= convert(T, 𝔸.B)
convert(::Type{InfinitesimalGenerator{U, V, C}}, M) where {U, V, C} = convert(BandedMatrix{U, V, C}, M)



#========================================================================================

Compute the operator
𝔸f = v_0 * f + v1 * ∂(f) + 0.5 * v2 * ∂∂(f)
𝔸'g = v_0 * g - ∂(v1 * g) + 0.5 * ∂∂(v2 * g)

========================================================================================#

function InfinitesimalGenerator(x::AbstractVector, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    InfinitesimalGenerator!(𝔸, make_Δ(x), v0, v1, v2)
end

function InfinitesimalGenerator!(𝔸, Δ, v0::AbstractVector, v1::AbstractVector, v2::AbstractVector)
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
    # Make sure each row sums to zero. Important in some cases: for isntance, otherwise cannot find sdf decomposition in GP model
    c = sum(𝔸, dims = 2)
    for i in 1:n
        𝔸[i, i] += v0[i] - c[i]
    end
    return InfinitesimalGenerator(𝔸)
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

Compute the principal eigenvector and eigenvalue of an operator

========================================================================================#
clean_eigenvalue(η::Union{Nothing, Real}) = η
function clean_eigenvalue(η::Complex)
    if abs(imag(η) .>= eps())
        @warn "Principal Eigenvalue has some imaginary part $(η)"
    end
    real(η)
end
clean_eigenvector_left(::Nothing) = nothing
clean_eigenvector_left(vl::Vector) = abs.(vl) ./ sum(abs.(vl))
clean_eigenvector_right(::Nothing) = nothing
clean_eigenvector_right(vr::Vector) = abs.(vr)



function principal_eigenvalue(T::InfinitesimalGenerator; method = :krylov, eigenvector = :right)
    η = nothing
    if method == :krylov
        vl, η, vr = principal_eigenvalue_krylov(T; eigenvector = eigenvector)
        if η == nothing
            @warn "Krylov Methods Failed"
        end
    end
    if η == nothing
        # use SuiteSparse maybe? LU decomposition sometimes?
        vl, η, vr = principal_eigenvalue_BLAS(convert(Matrix{Float64}, T); eigenvector = eigenvector)
    end
    return clean_eigenvector_left(vl), clean_eigenvalue(η), clean_eigenvector_right(vr)
end

# I could also use Arpack.eigs but it seems slower
function principal_eigenvalue_krylov(T::InfinitesimalGenerator; eigenvector = :right)
    vl, η, vr = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        vals, vecs, info = KrylovKit.eigsolve(T, 1, :LR, maxiter = size(T, 1))
        if info.converged > 0
            η = vals[1]
            vr = vecs[1]
        end
    end
    if eigenvector ∈ (:left, :both)
        vals, vecs, info = KrylovKit.eigsolve(T', 1, :LR, maxiter = size(T, 1))
        if info.converged > 0
            η = vals[1]
            vl = vecs[1]
        end
    end 
    return vl, η, vr
end

function principal_eigenvalue_BLAS(T::InfinitesimalGenerator; eigenvector = :right)
    vl, η, vr = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        e = eigen(T)
        _, out = findmax(real.(e.values))
        η = e.values[out]
        vr = e.vectors[:, out]
    end
    if eigenvector ∈ (:left, :both)
        e = eigen(copy(T'))
        _, out = findmax(real.(e.values))
        η = e.values[out]
        vl = e.vectors[:, out]
    end 
    return vl, η, vr
end

#========================================================================================

Feynman Kac

========================================================================================#

function feynman_kac_backward(𝔸::InfinitesimalGenerator; 
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


# Compute u(x, t)= E[∫0^t e^{-∫0^s V(x_τ)dτ}f(x_s)ds + e^{-∫0^tV(x_τ)dτ} ψ(x_t)|x_0 = x]
function feynman_kac_forward(𝔸::InfinitesimalGenerator; 
	t::AbstractVector = range(0, 100, step = 1/12), 
	ψ::AbstractVector = ones(size(𝔸, 1)), 
	f::AbstractVector = zeros(size(𝔸, 1)), 
	V::AbstractVector = zeros(size(𝔸, 1)))
    u = feynman_kac_backward(𝔸; ψ = ψ, t = - reverse(t), f = f, V = V)
    return u[:,end:-1:1]
end


#========================================================================================

Compute generator 𝔸f = E[df(x)]
where x is a diffusion process
dx = μx dt + σx dZ_t

========================================================================================#

function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)
    InfinitesimalGenerator(x, zeros(length(x)), μx, 0.5 * σx.^2)
end

# Stationary Distribution of x
function stationary_distribution(𝔸::InfinitesimalGenerator)
    principal_eigenvalue(𝔸; eigenvector = :left)[1]
end
function stationary_distribution(𝔸::InfinitesimalGenerator, δ, ψ)
    clean_eigenvector_left((δ * I - adjoint(𝔸)) \ (δ * ψ))
end


#========================================================================================

Compute generator 𝔸f = E[d(Mf(x))]
where x is a diffusive process
dx = μx dt + σx dZt
and M_t is a multiplicative functional
dMt/Mt = μM dt + σM dZt

========================================================================================#

function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector)
    InfinitesimalGenerator(x, μM, σM .* σx .+ μx, 0.5 * σx.^2)
end

# Compute Hansen Scheinkmann decomposition M = e^{ηt}f(x_t)W_t
function hansen_scheinkman(𝔸::InfinitesimalGenerator)
	principal_eigenvalue(𝔸; eigenvector = :right)[2:3]
end

##############################################################################
##
## Exported methods and types 
##
##############################################################################
export InfinitesimalGenerator,
principal_eigenvalue,
feynman_kac_backward,
feynman_kac_forward,
generator,
stationary_distribution,
hansen_scheinkman
end