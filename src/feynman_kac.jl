"""
    feynman_kac(𝕋 [; t, f, ψ, v]) 

With direction = :backward
Solve the following PDE:
u(x, t[end]) = ψ(x)
0 = u_t + 𝕋u - v(x, t)u + f(x, t)
Equivalently, in integral form, 
u(x, t) = E[∫_t^T e^{-∫_t^s v(x_u) du} f(x_s)ds + \int_t^t e^{-\int_t^T v(x_u)du} ψ(x_T)|x_t = x]
(notations are from the wikipedia article for Feynman–Kac formula)

With direction = :forward
Solve the following PDE:
u(x, t[1]) = ψ(x)
u_t = 𝕋u - v(x, t)u + f(x, t)
Equivalently, in integral form, 
u(x, t) = E[∫_0^t e^{-∫_0^s v(x_u) du} f(x_s)ds + \int_0^t e^{-\int_0^t v(x_u)du} ψ(x_t)|x_0 = x]

The function returns a matrix of size(length(f), length(t))
"""
function feynman_kac(𝕋; 
    t::AbstractVector = range(0, 100, step = 1/12), 
    f::Union{AbstractVector, AbstractMatrix} = zeros(size(𝕋, 1)), 
    ψ::AbstractVector = ones(size(𝕋, 1)),
    v::Union{AbstractVector, AbstractMatrix} = zeros(size(𝕋, 1)),
    direction= :backward)
    if direction == :backward
        u = zeros(size(𝕋, 1), length(t))
        u[:, end] = ψ
        if isa(f, AbstractVector) && isa(v, AbstractVector)
            if isa(t, AbstractRange)
                dt = step(t)
                B = factorize(I + (Diagonal(v) - 𝕋) * dt)
                for i in (length(t)-1):(-1):1
                    ψ = ldiv!(B, u[:, i+1] .+ f .* dt)
                    u[:, i] = ψ
                end
            else
                for i in (length(t)-1):(-1):1
                    dt = t[i+1] - t[i]
                    B = I + (Diagonal(v) - 𝕋) * dt
                    u[:, i] = B \ (u[:, i+1] .+ f .* dt)
                end
            end
        elseif isa(f, AbstractMatrix) && isa(v, AbstractMatrix)
            for i in (length(t)-1):(-1):1
                dt = t[i+1] - t[i]
                B = I + (Diagonal(view(v, :, i)) - 𝕋) * dt
                u[:, i] = B \ (u[:, i+1] .+ f[:, i] .* dt)
            end
        else
            error("f and v must be both AbstractVectors or both AbstractMatrices")
        end
        return u
    elseif direction == :forward
        u = feynman_kac(𝕋; t = - reverse(t), ψ = ψ, f = f, v = v, direction = :backward)
        return u[:,end:-1:1]
    else
        error("Direction must be :backward or :forward")
    end
end