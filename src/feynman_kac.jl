
#========================================================================================

Compute u(x_t, t)= E[∫t^T e^{-∫ts V(x_τ)dτ}f(x_s)ds + e^{-∫tTV(x_τ)dτ} ψ(x_T)|x_0 = x]

========================================================================================#

function feynman_kac_backward(x, μx, σx; ψ::AbstractVector, t::AbstractVector = range(0, 100, step = 1/12), f::AbstractVector = ones(length(x)), V::AbstractVector = zeros(length(x)))
    u = zeros(length(x), length(t))
    u[:, length(t)] = ψ
    ψ = deepcopy(ψ)
    A = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    if isa(t, StepRangeLen)
        dt = t[2] - t[1]
        A = factorize(I - build_operator!(A, x, V .* dt, μx .* dt, 0.5 .* σx.^2 .* dt)')
    end
    for i in (length(t)-1):(-1):1
        dt = t[i+1] - t[i]
        if !isa(t, StepRangeLen)
            A = I - build_operator!(A, x, V .* dt, μx .* dt, 0.5 .* σx.^2 .* dt)'
        end
        ψ = ldiv!(A, ψ .+ f .* dt)
        u[:, i] = ψ
    end
    return t, u
end

#========================================================================================

Compute u(x_t, T)= E[∫t^T e^{-∫ts V(x_τ)dτ}f(x_s)ds + e^{-∫tTV(x_τ)dτ} ψ(x_T)|x_t = x]

========================================================================================#

function feynman_kac_forward(x, μx, σx; ψ::AbstractVector, t::AbstractVector = range(0, 100, step = 1/12), f::AbstractVector = ones(length(x)), V::AbstractVector = zeros(length(x)))
    t, u feynman_kac_backward(x, μx, σx; ψ = ψ, t = t, f = f, V = V)
    reverse(t), X[:,end:-1:1]
end

#========================================================================================

Compute u(x_t, t) = E[∫t^T e^{-∫ts V(x_τ, τ)dτ}f(x_s, s)ds + e^{-∫tT V(x_τ, τ)dτ}ψ(x_T)|x_t = x]

========================================================================================#

function feynman_kac_backward(x, μx, σx; ψ::AbstractVector, t::AbstractVector = range(0, 100, step = 1/12), f::AbstractMatrix = ones(length(x), length(t)), V::AbstractVector = zeros(length(x), length(t)))
    u = zeros(length(x), length(t))
    u[:, length(t)] = ψ
    ψ = deepcopy(ψ)
    A = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    for i in (length(t)-1):(-1):1
        dt = t[i+1] - t[i]
        A = I - build_operator!(A, x, V[:, t] .* dt, μx .* dt, 0.5 .* σx.^2 .* dt)'
        ψ = ldiv!(A, ψ .+ f[:, i] .* dt)
        u[:, i] = ψ
    end
    return t, u
end
