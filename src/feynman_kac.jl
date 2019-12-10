"""
With direction = :backward
Solve the PDE backward in time
u(x, T) = ψ(x)
0 = u_t + 𝔸u_t - V(x, t)u +  f(x, t)


With direction = :forward
Solve the PDE forward in time
u(x, 0) = ψ(x)
u_t = 𝔸u - V(x)u + f(x)
"""
function feynman_kac(𝔸::AbstractMatrix; 
    t::AbstractVector = range(0, 100, step = 1/12), 
    ψ::AbstractVector = ones(size(𝔸, 1)), 
    f::Union{AbstractVector, AbstractMatrix} = zeros(size(𝔸, 1)), 
    V::Union{AbstractVector, AbstractMatrix} = zeros(size(𝔸, 1)),
    direction= :backward)
    if direction == :backward
        u = zeros(size(𝔸, 1), length(t))
        u[:, end] = ψ
        if isa(f, AbstractVector) && isa(V, AbstractVector)
            if isa(t, AbstractRange)
                dt = step(t)
                𝔹 = factorize(I + (Diagonal(V) - 𝔸) * dt)
                for i in (length(t)-1):(-1):1
                    ψ = ldiv!(𝔹, u[:, i+1] .+ f .* dt)
                    u[:, i] = ψ
                end
            else
                for i in (length(t)-1):(-1):1
                    dt = t[i+1] - t[i]
                    𝔹 = I + (Diagonal(V) - 𝔸) * dt
                    u[:, i] = 𝔹 \ (u[:, i+1] .+ f .* dt)
                end
            end
        elseif isa(f, AbstractMatrix) && isa(V, AbstractMatrix)
            for i in (length(t)-1):(-1):1
                dt = t[i+1] - t[i]
                𝔹 = I + (Diagonal(view(V, :, i)) - 𝔸) * dt
                u[:, i] = 𝔹 \ (u[:, i+1] .+ f[:, i] .* dt)
            end
        else
            error("f and V must be Vectors or Matrices")
        end
        return u
    elseif direction == :forward
        u = feynman_kac(𝔸; t = - reverse(t), ψ = ψ, f = f, V = V, direction = :backward)
        return u[:,end:-1:1]
    else
        error("Direction must be :backward or :forward")
    end
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



""" 
If direction = :forward
compute `E[M_t ψ(x_t)|x_0 = x]`
"""
function feynman_kac(M::MultiplicativeFunctional; kwargs...)
    feynman_kac(generator!(M); kwargs...)
end

