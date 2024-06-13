"""
    feynman_kac(𝕋, ts; f =  zeros(size(𝕋, 1)), ψ =  zeros(size(𝕋, 1)), v = zeros(size(𝕋, 1)), direction = :backward]) 

𝕋 should be a matrix 
ts should be a grid of time on which to solve the PDE

With direction = :backward, returns the solution of the PDE:
u(x, t[end]) = ψ(x)
0 = u_t + 𝕋u - v(x, t)u + f(x, t)
Or, equivalently, in integral form, 
u(x, t) = E[∫_t^T e^{-∫_t^s v(x_u) du} f(x_s)ds + ∫_t^t e^{-∫_t^T v(x_u)du} ψ(x_T)|x_t = x]
(notations are from the wikipedia article for Feynman–Kac formula)

With direction = :forward, , returns the solution of the PDE:
u(x, t[1]) = ψ(x)
u_t = 𝕋u - v(x, t)u + f(x, t)
Or, equivalently, in integral form, 
u(x, t) = E[∫_0^t e^{-∫_0^s v(x_u) du} f(x_s)ds + ∫_0^t e^{-∫_0^t v(x_u)du} ψ(x_t)|x_0 = x]

The PDE is solved using Euler method with implicit time steps
"""
function feynman_kac(𝕋, ts; 
    f::Union{AbstractVector, AbstractMatrix} = zeros(size(𝕋, 1)), 
    ψ::AbstractVector = zeros(size(𝕋, 1)),
    v::Union{AbstractVector, AbstractMatrix} = zeros(size(𝕋, 1)),
    direction= :backward)
    size(𝕋, 1) == size(𝕋, 2) || throw(DimensionMismatch(), "𝕋 must be square matrix")
    size(𝕋, 1) == size(f, 1) || throw(DimensionMismatch(), "𝕋 and f should have the same number of rows")
    size(𝕋, 1) == length(ψ) || throw(DimensionMismatch(), "𝕋 and ψ should have the same number of rows")
    size(𝕋, 1) == size(v, 1) || throw(DimensionMismatch(), "𝕋 and v should have the same number of rows")
    size(f, 2) ∈ (1, length(ts)) ||  throw(DimensionMismatch(), "The number of columns in f should equal the length of ts")
    size(v, 2) ∈ (1, length(ts)) ||  throw(DimensionMismatch(), "The number of columns in f should equal the length of ts")
    direction ∈ (:forward, :backward) || throw(ArgumentError(), "Direction must be :backward or :forward")
    if ndims(f) == 2 && ndims(v) == 1
        v = hcat([v for _ in 1:size(f, 2)])
    elseif ndims(f) == 1 && ndims(v) == 2
        f = hcat([f for _ in 1:size(v, 2)])
    end
    if direction == :forward
        # direction is forward
        u = feynman_kac(𝕋, - reverse(ts); ψ = ψ, f = f, v = v, direction = :backward)
        return u[:,end:-1:1]
    else
        # direction is backward
        u = zeros(size(𝕋, 1), length(ts))
        u[:, end] = ψ
        if ndims(f) == 1
            # f and v are vectors
            if isa(ts, AbstractRange)
                # constant time step
                dt = step(ts)
                B = factorize(I + (Diagonal(v) - 𝕋) * dt)
                for i in (length(ts)-1):(-1):1
                    ψ = ldiv!(B, u[:, i+1] .+ f .* dt)
                    u[:, i] = ψ
                end
            else
                # non-constant time step
                for i in (length(ts)-1):(-1):1
                    dt = ts[i+1] - ts[i]
                    B = I + (Diagonal(v) - 𝕋) * dt
                    u[:, i] = B \ (u[:, i+1] .+ f .* dt)
                end
            end
        else
            # f and v are matrices
            for i in (length(ts)-1):(-1):1
                dt = ts[i+1] - ts[i]
                B = I + (Diagonal(view(v, :, i)) - 𝕋) * dt
                u[:, i] = B \ (u[:, i+1] .+ f[:, i] .* dt)
            end
        end
        return u
    end
end