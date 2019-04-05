#========================================================================================

Compute generator 𝔸f = E[d(Mf(x))]
where x is a diffusive process
dx = μx dt + σx dZt
and M_t is a multiplicative functional
dMt/Mt = μM dt + σM dZt

========================================================================================#

function generator(x::AbstractVector, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector)
    𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    Δ = make_Δ(x)
    generator!(𝔸, Δ, μx, σx, μM, σM)
end

function generator!(𝔸::AbstractMatrix, Δ, μx::AbstractVector, σx::AbstractVector, μM::AbstractVector, σM::AbstractVector)
    operator!(𝔸, Δ, μM, σM .* σx .+ μx, 0.5 * σx.^2)
end

#========================================================================================

Compute Hansen Scheinkmann decomposition M = e^{ηt}f(x_t)W_t
where x is a diffusive process
dx = μx dt + σx dZt
and M_t is a multiplicative functional
dMt/Mt = μM dt + σM dZt

The function returns g, η, f

========================================================================================#
function hansen_scheinkman(x, μx, σx, μM, σM; method = :krylov, eigenvector = :right)
    𝔸 = BandedMatrix(Zeros(length(x), length(x)), (1, 1))
    Δ = make_Δ(x)
    hansen_scheinkman!(𝔸, Δ, μx, σx, μM, σM; method = method, eigenvector = eigenvector)
end

function hansen_scheinkman!(𝔸, Δ, μx, σx, μM, σM; method = :krylov, eigenvector = :right)
    generator!(𝔸, Δ, μx, σx, μM, σM)
    principal_eigenvalue(𝔸; method = method, eigenvector = eigenvector)
end

#========================================================================================

Compute u(x, t) = E[M_tψ(x_t)|x_0 = x] using Implicit Feynman Kac
where x is a diffusive process
dx = μx dt + σx dZt
and M_t is a multiplicative functional
dMt/Mt = μM dt + σM dZt

========================================================================================#

function feynman_kac_forward(x, μx, σx, μM, σM; t::AbstractVector = range(0, 100, step = 1/12), ψ = ones(length(x)))
    feynman_kac_forward(x, μx .+ σM .* σx, σx; t = t, ψ = ψ, V = - μM)
end

#========================================================================================

Compute ϵ(x, T) = σD(x) * (σM + σE[M_T | X_0 = x])
where x is a diffusive process
dx = μx dt + σx dZt
and M_t is a multiplicative functional
dMt/Mt = μM dt + σM dZt

========================================================================================#

function impulse_response(x, μx, σx, μM, σM; t::AbstractVector = range(0, 100, step = 1/12),  σD = ones(length(x)))
    u = feynman_kac_forward(x, μx, σx, μM, σM; t = t)
    for i in 1:length(t)
        u[:, i] = σD .* (σM .+ _derive(u[:, i], x, μx) ./ u[:, i] .* σx)
    end
    return u
end

function _derive(f::AbstractVector, x::AbstractVector, μx::AbstractVector)
    out = similar(f)
    n = length(f)
    for i in 1:n
        if μx[i] >= 0
            out[i] = (f[min(i+1, n)] - f[i]) / (x[min(i+1, n)] - x[i])
        else
            out[i] = (f[i] - f[max(i-1, 1)]) / (x[i] - x[max(i-1, 1)])
        end
    end
    return out
end

