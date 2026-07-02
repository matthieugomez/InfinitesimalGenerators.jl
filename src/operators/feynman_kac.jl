"""
    feynman_kac(𝕋, ts; f =  zeros(size(𝕋, 1)), ψ =  zeros(size(𝕋, 1)), v = zeros(size(𝕋, 1)), direction = :backward)

𝕋 should be a matrix
ts should be a grid of time on which to solve the PDE

With direction = :backward, returns the solution of the PDE:
u(x, t[end]) = ψ(x)
0 = u_t + 𝕋u - v(x, t)u + f(x, t)
Or, equivalently, in integral form,
u(x, t) = E[∫_t^T e^{-∫_t^s v(x_u) du} f(x_s)ds + e^{-∫_t^T v(x_u)du} ψ(x_T)|x_t = x]
(notations are from the wikipedia article for Feynman–Kac formula)

With direction = :forward, returns the solution of the PDE:
u(x, t[1]) = ψ(x)
u_t = 𝕋u - v(x, t)u + f(x, t)
Or, equivalently, in integral form,
u(x, t) = E[∫_0^t e^{-∫_0^s v(x_u) du} f(x_s)ds + e^{-∫_0^t v(x_u)du} ψ(x_t)|x_0 = x]

The PDE is solved using Euler method with implicit time steps
"""
function _flatten_state_vector_argument(X::MarkovProcess, x, label::Symbol)
    length(x) == length(X) ||
        throw(DimensionMismatch("`$label` has length $(length(x)) but process has length $(length(X))"))
    return vec(x)
end

function _flatten_state_time_argument(X::MarkovProcess, x, ts, label::Symbol)
    state_shape = size(X)
    if x isa AbstractVector
        length(x) == length(X) ||
            throw(DimensionMismatch("`$label` has length $(length(x)) but process has length $(length(X))"))
        return x
    elseif size(x) == state_shape
        return vec(x)
    elseif x isa AbstractMatrix && size(x, 1) == length(X) && size(x, 2) in (1, length(ts))
        return x
    elseif length(size(x)) == length(state_shape) + 1 &&
            size(x)[1:end - 1] == state_shape &&
            size(x, ndims(x)) in (1, length(ts))
        return reshape(x, length(X), size(x, ndims(x)))
    else
        throw(DimensionMismatch(
            "`$label` must have shape $(state_shape), length $(length(X)), " *
            "or time-varying shape ($(state_shape)..., $(length(ts)))"))
    end
end

function feynman_kac(X::MarkovProcess, ts; f = nothing, ψ = nothing, v = nothing, direction = :backward)
    𝕋 = generator(X)
    f_flat = f === nothing ? zeros(eltype(𝕋), length(X)) :
        _flatten_state_time_argument(X, f, ts, :f)
    ψ_flat = ψ === nothing ? zeros(eltype(𝕋), length(X)) :
        _flatten_state_vector_argument(X, ψ, :ψ)
    v_flat = v === nothing ? zeros(eltype(𝕋), length(X)) :
        _flatten_state_time_argument(X, v, ts, :v)
    u = feynman_kac(𝕋, ts; f = f_flat, ψ = ψ_flat, v = v_flat, direction = direction)
    return _reshape_state_time_output(X, u)
end

function feynman_kac(𝕋, ts;
    f::Union{AbstractVector, AbstractMatrix} = zeros(eltype(𝕋), size(𝕋, 1)),
    ψ::AbstractVector = zeros(eltype(𝕋), size(𝕋, 1)),
    v::Union{AbstractVector, AbstractMatrix} = zeros(eltype(𝕋), size(𝕋, 1)),
    direction= :backward)
    size(𝕋, 1) == size(𝕋, 2) || throw(DimensionMismatch("𝕋 must be square matrix"))
    size(𝕋, 1) == size(f, 1) || throw(DimensionMismatch("𝕋 and f should have the same number of rows"))
    size(𝕋, 1) == length(ψ) || throw(DimensionMismatch("𝕋 and ψ should have the same number of rows"))
    size(𝕋, 1) == size(v, 1) || throw(DimensionMismatch("𝕋 and v should have the same number of rows"))
    size(f, 2) ∈ (1, length(ts)) ||  throw(DimensionMismatch("The number of columns in f should equal the length of ts"))
    size(v, 2) ∈ (1, length(ts)) ||  throw(DimensionMismatch("The number of columns in v should equal the length of ts"))
    direction ∈ (:forward, :backward) || throw(ArgumentError("Direction must be :backward or :forward"))
    if ndims(f) == 2 && ndims(v) == 1
        v = repeat(v, 1, size(f, 2))
    elseif ndims(f) == 1 && ndims(v) == 2
        f = repeat(f, 1, size(v, 2))
    end
    if direction == :forward
        # direction is forward
        f_reverse = ndims(f) == 2 ? @view(f[:, end:-1:1]) : f
        v_reverse = ndims(v) == 2 ? @view(v[:, end:-1:1]) : v
        u = feynman_kac(𝕋, - reverse(ts); ψ = ψ, f = f_reverse, v = v_reverse, direction = :backward)
        return u[:,end:-1:1]
    else
        # direction is backward
        T = float(promote_type(eltype(𝕋), eltype(f), eltype(ψ), eltype(v), eltype(ts)))
        u = zeros(T, size(𝕋, 1), length(ts))
        u[:, end] = ψ
        if ndims(f) == 1
            # f and v are vectors
            if isa(ts, AbstractRange)
                # constant time step
                dt = step(ts)
                B = factorize(I + (Diagonal(v) - 𝕋) * dt)
                for i in (length(ts)-1):(-1):1
                    u[:, i] = B \ (u[:, i+1] .+ f .* dt)
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
