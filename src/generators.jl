#========================================================================================

For a Markov Process x:
dx = μx dt + σx dZ_t

========================================================================================#

mutable struct MarkovProcess
    x::AbstractVector{<:Real}
    μx::AbstractVector{<:Real}
    σx::AbstractVector{<:Real}
    𝔸::Tridiagonal
    Δ::Tuple{<:AbstractVector, <:AbstractVector, <:AbstractVector, <:AbstractVector}
end

function MarkovProcess(x::AbstractVector{<:Real}, μx::AbstractVector{<:Real}, σx::AbstractVector{<:Real})
    length(x) == length(μx) || error("Vector for grid, drift, and volatility should have the same size")
    length(μx) == length(σx) || error("Vector for grid, drift, and volatility should have the same size")
    n = length(x)
    𝔸 = Tridiagonal(zeros(n-1), zeros(n), zeros(n-1))
    Δ = make_Δ(x)
    MarkovProcess(x, μx, σx, 𝔸, Δ)
end

Base.length(x::MarkovProcess) = length(x.x)

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

# it's important to take 1e-6 to have the right tail index of multiplicative functional (see tests)
function OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, 
    xmin = quantile(Normal(xbar, σ / sqrt(2 * κ)), p), xmax = quantile(Normal(xbar, σ / sqrt(2 * κ)), 1 - p))
    x = range(xmin, stop = xmax, length = length)
    μx = κ .* (xbar .- x)
    σx = σ .* Ones(Base.length(x))
    MarkovProcess(x, μx, σx)
end

function CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0, p = 1e-10, length = 100, α = 2 * κ * xbar / σ^2, β = σ^2 / (2 * κ), xmin = quantile(Gamma(α, β), p), xmax = quantile(Gamma(α, β), 1 - p), pow = 2)
    x = range(xmin^(1/pow), stop = xmax^(1/pow), length = length).^pow
    μx = κ .* (xbar .- x)
    σx = σ .* sqrt.(x)
    MarkovProcess(x, μx, σx)
end

# Compute generator 𝔸f = E[df(x)]
function generator!(x::MarkovProcess)
    operator!(x.𝔸, x.Δ, Zeros(length(x.x)), x.μx, 0.5 * x.σx.^2)
end

function generator(x::MarkovProcess)
    deepcopy(generator!(x))
end

# Stationary Distribution of x
function stationary_distribution(x::MarkovProcess)
    g, η, _ = principal_eigenvalue(generator!(x); which = :SM, eigenvector = :left)
    abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
    return g
end

# Stationary Distribution of x with death rate δ and reinjection ψ
function stationary_distribution(x::MarkovProcess, δ::Number, ψ::AbstractVector{<:Number})
    clean_eigenvector_left((δ * I - generator!(x)') \ (δ * ψ))
end

function ∂(x::MarkovProcess, f::AbstractVector)
	operator!(x.𝔸, x.Δ, Zeros(length(x.x)), Ones(length(x.x)), Zeros(length(x.x))) * f
end

#========================================================================================

For a multiplicative functional M:
dM/M = μM(x) dt + σM(x) dZt

========================================================================================#

mutable struct MultiplicativeFunctional
    x::MarkovProcess
    μM::AbstractVector{<:Number}
    σM::AbstractVector{<:Number}
    ρ::Number
    δ::Number
end

function MultiplicativeFunctional(x::MarkovProcess, μM::AbstractVector{<:Number}, σM::AbstractVector{<:Number}; ρ::Number = 0.0, δ::Number = 0.0)
    length(x.x) == length(μM) || error("Vector for grid and μM should have the same size")
    length(x.x) == length(σM) || error("Vector for grid and σM should have the same size")
    MultiplicativeFunctional(x, μM, σM, ρ, δ)
end

MarkovProcess(M::MultiplicativeFunctional) = M.x
Base.length(M::MultiplicativeFunctional) = length(M.x)
# Generator for long run CGF
function generator!(M::MultiplicativeFunctional, ξ = 1.0)
    operator!(M.x.𝔸, M.x.Δ, ξ .* M.μM .+ 0.5 * ξ * (ξ - 1) .* M.σM.^2 .- M.δ,  M.x.μx .+ ξ .* M.σM .* M.ρ .* M.x.σx, 0.5 * M.x.σx.^2)
end
function generator(M::MultiplicativeFunctional, ξ = 1.0)
    deepcopy(generator!(M, ξ))
end

# ξ -> lim log(E[M^\xi]) / t
function cgf_longrun(M::MultiplicativeFunctional; which = :LR, eigenvector = :right, r0 = Ones(length(M.x)))
    ξ -> principal_eigenvalue(generator!(M, ξ), which = which, eigenvector = eigenvector, r0 = r0)
end

# Compute Hansen Scheinkmann decomposition M_t= e^{ηt}f(x_t)\hat{M}_t
function hansen_scheinkman(M::MultiplicativeFunctional; which = :LR, eigenvector = :right)
    cgf_longrun(M, eigenvector = eigenvector)(1.0)
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
function tail_index(M::MultiplicativeFunctional; which = :SM, xatol = 1e-2, verbose = false, kwargs...)
    out = 0.0
    r0 = ones(length(M.x))
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