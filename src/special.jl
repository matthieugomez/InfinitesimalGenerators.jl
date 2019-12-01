struct MarkovProcess
    x::AbstractVector{<:Real}
    μx::AbstractVector{<:Real}
    σx::AbstractVector{<:Real}
end



function OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.1, σ = 1.0, length = 100)
    x = range(xbar - 5 * σ , xbar + 5 * σ, length = 100)
    μx .= κ .* (xbar .- x)
    σx .= σ .* ones(length(x))
    (x, μx, σx)
end


function CoxIngersollRoss(; xbar = 0.1, κ = 0.1, σ = 1.0, length = 100)
        x = range(μ - 5 * σ , μ + 5 * σ, length = 100)
        α = 2 * κ * xbar / σ^2
        β = σ^2 / (2 * κ)
        vmin = quantile(Gamma(α, β), 0.025)
        vmax = quantile(Gamma(α, β), 0.975)
        μx .= κ .* (μ .- x)
        σx .= σ .* sqrt.(x)
        (x, μx, σx)
end
