
#========================================================================================

Stationary Distribution with one state variable

========================================================================================#
#now there are still two issues
#1. Does not satisfy walras law. Or mathematically does not satisfy IPP ∑ μ.g = ∑ a.Ag. 
# 1.1. First part due to drift if not positive at left boundary or not negative ar right boundary In the case drift is positive, there is a remaning term μ_NdG(a_N) To satisfy it, do amax super super high (intuitively, x high enough so that cutting behavior at the top does not matter for aggregate as g(x)x -> 0)
#1.2 Second part is due to volatility. Note that it requires to put invΔx[i] for central derivative, which is different with the formula in Moll notes
#2. A g can be negative when updating forward. Use implicit scheme


function stationary_distribution(grid, μ::Vector{T}, σ2::Vector{T}) where {T <: Number}
    A = build_operator(grid, zero(grid), μ, 0.5 * σ2)
    _, η, density = principal_eigenvalue(A; eigenvector = :right)
    abs.(density ./ sum(abs.(density)))
end


function stationary_distribution(grid, μ::Vector{T}, σ2::Vector{T}, δ, ψ) where {T <: Number}
    A = build_operator(grid, zero(grid), μ, 0.5 * σ2)
    density = (δ * I - A) \ (δ * ψ)
    density .= abs.(density ./ sum(abs.(density)))
end