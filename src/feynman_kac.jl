
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
                ψ = 𝔹 \ (u[:, i+1] .+ f .* dt)
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
    t::AbstractVector = range(0, 100, step = 1/12), kwargs...)
    u = feynman_kac_backward(𝔸; t = - reverse(t), kwargs...)
    return u[:,end:-1:1]
end
