
abstract type MarkovProcess end

# This type sould define generator(), which returns a transition matrix 𝕋 such that
# 𝕋f = lim_{t→0} E[f(x_t)|x_0=x]/t


abstract type AbstractGenerator <: AbstractMatrix end
struct InfinitesimalGenerator{T} <: AbstractGenerator
    M::T
    function InfinitesimalGenerator(M)
        all(x ≈ 0 for y in sum(M, dims = 2)) || throw(ArgumentError("Rows do not sum up to one"))
        check_offdiagonal_negative(M)
        return InfinitesimalGenerator(M)
    end
end

function check_offdiagonal_negative(A::AbstractMatrix)
    @inbounds for ι in CartesianIndices(A)
        if ι[1] ≠ ι[2]
            a = A[ι]
            A[ι] ≥ 0 || throw(ArgumentError("Element at $(ι) is negative."))
        end
    end
end

function check_offdiagonal_negative(A::SparseMatrixCSC)
    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)
    for j in 1:n
        for k in nzrange(A, j)
            i = rows[k]
            if i ≠ j
                vals[k] || throw(DomainError(a, "Element at $(CartesianIndex(i,j)) is negative."))
            end
        end
    end
    d
end

function check_offdiagonal_negative(A::Tridiagonal)
    all(x <= 0 for x in A.dl) || throw(ArgumentError("off-diagonal elements are positive"))
    all(x <= 0 for x in A.du) || throw(ArgumentError("off-diagonal elements are positive"))
    d
end



for fun in (:conj, :copy, :real, :imag, :adjoint, :ishermitian, :issymmetric, :iszero, :isone, :istriu, :isdiag)
    @eval function ($func)(M::AbstractGenerator)
        $func(M.M)
    end
end

for fun in (:getindex, :setindex!, Base.replace_in_print_matrix)
        @eval function ($func)(M::AbstractGenerator, ...)
            $func(M.M)
        end
    end
end



"""
 computes the stationary distribution corresponding to the MarkovProcess X
"""
function stationary_distribution(X::MarkovProcess; δ = 0.0, ψ = Zeros(length(X.x)))
    δ >= 0 ||  throw(ArgumentError("δ needs to be positive"))
    if δ > 0
        g = abs.((δ * I - generator(X)') \ (δ * ψ))
    else
        η, g = principal_eigenvalue(generator(X)')
        abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
    end
    g ./ sum(g)
end


"""
    Returns the Diffusion Process `x_t` with SDE
    
        dx_t = μ(x_t) dt + σ(x_t) dZ_t

"""
mutable struct DiffusionProcess <: MarkovProcess
    x::AbstractVector{<:Real}
    μx::AbstractVector{<:Real}
    σx::AbstractVector{<:Real}
    function DiffusionProcess(x::AbstractVector{<:Real}, μx::AbstractVector{<:Real}, σx::AbstractVector{<:Real})
        length(x) == length(μx) == length(σx) || throw(ArgumentError("Vector for grid, drift, and volatility should have the same size"))
        new(x, μx, σx)
    end
end

state_space(X::DiffusionProcess) = X.x

"""
    Returns the discretized version of the infinitesimal generator of the Diffusion Process
    
        𝕋: f ⭌ lim 1/t * E[f(x_t)|x_0=x]
                 = μx * ∂f + 0.5 * σx^2 * ∂^2f

    defined on the set of functions f such that 
        
        ∂f(x) = 0 

    at the border of the state space

    The transpose of this operator corresponds to
        
        𝕋': g ⭌ v * g - ∂(μx * g) + 0.5 * ∂^2(σx^2 * g)

    defined on the set of functions g such that  
        
        -μx * g(x) + 0.5 * ∂(σx^2 * g) = 0

    at the border of state space
"""
generator(X::DiffusionProcess) = generator(X.x, X.μx, X.σx)

function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)
    n = length(x)
    𝕋 = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    @inbounds for i in 1:n
        Δxp = x[min(i, n-1)+1] - x[min(i, n-1)]
        Δxm = x[max(i-1, 1) + 1] - x[max(i-1, 1)]
        Δx = (Δxm + Δxp) / 2
        # upwinding to ensure off diagonals are posititive
        if μx[i] >= 0
            𝕋[i, min(i + 1, n)] += μx[i] / Δxp
            𝕋[i, i] -= μx[i] / Δxp
        else
            𝕋[i, i] += μx[i] / Δxm
            𝕋[i, max(i - 1, 1)] -= μx[i] / Δxm
        end
        𝕋[i, max(i - 1, 1)] += 0.5 * σx[i]^2 / (Δxm * Δx)
        𝕋[i, i] -= 0.5 * σx[i]^2 * 2 / (Δxm * Δxp)
        𝕋[i, min(i + 1, n)] += 0.5 * σx[i]^2 / (Δxp * Δx)
    end
    # ensure rows sum to zero with machine precision
    c = sum(𝕋, dims = 2)
    for i in 1:n
        𝕋[i, i] -= c[i]
    end
    return 𝕋
end


"""
    Returns the discretized version of the operator ∂
    
        δ: f ⭌ ∂f

"""
function ∂(X::DiffusionProcess)
    Diagonal(X.μx) \ generator(X.x, X.μx, Zeros(length(X.x)))
end


"""
    Returns the Ornstein Uhlenbeck process defined by the SDE
        
        dx_t = -κ * (x_t - xbar) * dt + σ * dZ_t

"""
function OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, 
    xmin = quantile(Normal(xbar, σ / sqrt(2 * κ)), p), xmax = quantile(Normal(xbar, σ / sqrt(2 * κ)), 1 - p))
    # it's important to take low p to have the right tail index of Additive functional
    if xmin > 0
        x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    else
        x = range(xmin, stop = xmax, length = length)
    end
    DiffusionProcess(x, κ .* (xbar .- x), σ * Ones(Base.length(x)))
end

"""
    Returns the Cox Ingersoll Ross process defined by the SDE
        
        dx_t = -κ * (x - xbar) * dt + σ * sqrt(x) * dZ_t

"""
function CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, α = 2 * κ * xbar / σ^2, β = σ^2 / (2 * κ), xmin = quantile(Gamma(α, β), p), xmax = quantile(Gamma(α, β), 1 - p), pow = 2)
    # check 0 is not attainable
    @assert (2 * κ * xbar) / σ^2 > 1
    x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    DiffusionProcess(x, κ .* (xbar .- x), σ .* sqrt.(x))
end

