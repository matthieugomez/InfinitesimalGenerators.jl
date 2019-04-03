
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
clean_density(::Nothing) = nothing
clean_density(v::Vector) = abs.(v) ./ sum(abs.(v))

clean_f(v::Vector) = abs.(v)
clean_f(::Nothing) = nothing

function principal_eigenvalue(T; method = :krylov, eigenvector = :right)
    η = nothing
    if method == :krylov
        vl, η, vr = principal_eigenvalue_krylov(T; eigenvector = eigenvector)
        if η == nothing
            @warn "Krylov Methods Failed"
        end
    end
    if η == nothing
        vl, η, vr = principal_eigenvalue_BLAS(convert(Matrix{Float64}, T); eigenvector = eigenvector)
    end
    return clean_f(vl), clean_eigenvalue(η), clean_density(vr)
end

function principal_eigenvalue_krylov(T; eigenvector = :right)
    vl, η, vr = nothing, nothing, nothing
    if eigenvector ∈ (:right, :both)
        vals, vecs, info = KrylovKit.eigsolve(T, 1, :LR)
        if info.converged > 0
            η = vals[1]
            vr = vecs[1]
        end
    end
    if eigenvector ∈ (:left, :both)
        vals, vecs, info = KrylovKit.eigsolve(T', 1, :LR)
        if info.converged > 0
            η = vals[1]
            vl = vecs[1]
        end
    end 
    return vl, η, vr
end


function principal_eigenvalue_BLAS(T; eigenvector = :right)
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

Stationary Distribution with one state variable

========================================================================================#
#now there are still two issues
#1. Does not satisfy walras law. Or mathematically does not satisfy IPP ∑ μ.g = ∑ a.Ag. 
# 1.1. First part due to drift if not positive at left boundary or not negative ar right boundary In the case drift is positive, there is a remaning term μ_NdG(a_N) To satisfy it, do amax super super high (intuitively, x high enough so that cutting behavior at the top does not matter for aggregate as g(x)x -> 0)
#1.2 Second part is due to volatility. Note that it requires to put invΔx[i] for central derivative, which is different with the formula in Moll notes
#2. A g can be negative when updating forward. Use implicit scheme


function stationary_distribution(x::AbstractVector, μ::AbstractVector, σ::AbstractVector)
    A = build_operator(x, zero(x), μ, 0.5 * σ.^2)
    _, _, density = principal_eigenvalue(A; eigenvector = :right)
    clean_density(density)
end


function stationary_distribution(x::AbstractVector, μ::AbstractVector, σ::AbstractVector, δ, ψ)
    A = build_operator(x, zero(x), μ, 0.5 * σ.^2)
    density = (δ * I - A) \ (δ * ψ)
    clean_density(density)
end