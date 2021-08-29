"""
With direction = :backward
Solve the PDE backward in time
u(x, t[end]) = ψ(x)
0 = u_t + Tu_t - V(x, t)u +  f(x, t)


With direction = :forward
Solve the PDE forward in time
u(x, t[1]) = ψ(x)
u_t = Tu - V(x)u + f(x)
"""
function feynman_kac(𝕋; 
    t::AbstractVector = range(0, 100, step = 1/12), 
    f::Union{AbstractVector, AbstractMatrix} = zeros(size(𝕋, 1)), 
    ψ::AbstractVector = ones(size(𝕋, 1)),
    V::Union{AbstractVector, AbstractMatrix} = zeros(size(𝕋, 1)),
    direction= :backward)
    if direction == :backward
        u = zeros(size(𝕋, 1), length(t))
        u[:, end] = ψ
        if isa(f, AbstractVector) && isa(V, AbstractVector)
            if isa(t, AbstractRange)
                dt = step(t)
                B = factorize(I + (Diagonal(V) - 𝕋) * dt)
                for i in (length(t)-1):(-1):1
                    ψ = ldiv!(B, u[:, i+1] .+ f .* dt)
                    u[:, i] = ψ
                end
            else
                for i in (length(t)-1):(-1):1
                    dt = t[i+1] - t[i]
                    B = I + (Diagonal(V) - 𝕋) * dt
                    u[:, i] = B \ (u[:, i+1] .+ f .* dt)
                end
            end
        elseif isa(f, AbstractMatrix) && isa(V, AbstractMatrix)
            for i in (length(t)-1):(-1):1
                dt = t[i+1] - t[i]
                B = I + (Diagonal(view(V, :, i)) - 𝕋) * dt
                u[:, i] = B \ (u[:, i+1] .+ f[:, i] .* dt)
            end
        else
            error("f and V must be Vectors or Matrices")
        end
        return u
    elseif direction == :forward
        u = feynman_kac(𝕋; t = - reverse(t), ψ = ψ, f = f, V = V, direction = :backward)
        return u[:,end:-1:1]
    else
        error("Direction must be :backward or :forward")
    end
end