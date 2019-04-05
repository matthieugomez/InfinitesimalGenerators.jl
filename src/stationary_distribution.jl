

#========================================================================================

Stationary Distribution with one state variable

========================================================================================#
#now there are still two issues
#1. Does not satisfy walras law. Or mathematically does not satisfy IPP ∑ μ.g = ∑ a.Ag. 
# 1.1. First part due to drift if not positive at left boundary or not negative ar right boundary In the case drift is positive, there is a remaning term μ_NdG(a_N) To satisfy it, do amax super super high (intuitively, x high enough so that cutting behavior at the top does not matter for aggregate as g(x)x -> 0)
#1.2 Second part is due to volatility. Note that it requires to put invΔx[i] for central derivative, which is different with the formula in Moll notes
#2. A g can be negative when updating forward. Use implicit scheme

function compute_generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)
    𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    Δ = make_Δ(x)
    compute_generator!(𝔸, Δ, μx, σx)
end

function compute_generator!(𝔸::AbstractMatrix, Δ, μx::AbstractVector, σx::AbstractVector)
    build_operator!(𝔸, Δ, zeros(length(x)), μx, 0.5 * σx.^2)
end



function stationary_distribution(x::AbstractVector, μx::AbstractVector, σx::AbstractVector)
    𝔸 = compute_generator(x, μx, σx)
    density, _, _ = principal_eigenvalue(𝔸; eigenvector = :right)
end


function stationary_distribution(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, δ, ψ)
    𝔸 = compute_generator(x, μx, σx)
    density = (δ * I - 𝔸') \ (δ * ψ)
    clean_density(density)
end