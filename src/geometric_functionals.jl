#========================================================================================

Compute Hansen Scheinkmann decomposition M = e^{ηt}f(x_t)W_t
Return g, η, f

========================================================================================#
function compute_η(x, μx, σx, μM, σM; method = :krylov, eigenvector = :right)
    n = length(x)
    𝔸 = zeros(n, n)
    Δ = EconPDEs.make_Δ(x)
    compute_η!(𝔸, Δ, μx, σx, μM, σM; method = method, eigenvector = eigenvector)
end

function compute_η!(𝔸, Δ, μx, σx, μM, σM; method = :krylov, eigenvector = :right)
    build_operator!(𝔸, Δ, μM, σM .* σx .+ μx, 0.5 .* σx.^2)
    principal_eigenvalue(𝔸; method = method, eigenvector = eigenvector)
end


#========================================================================================

Compute u(x, T) = E[M_Tψ(x_T)|x_t = x] using Implicit Feynman Kac
where
dx = μx dt + σx dZ_t
and M_t is a geometric functional
dMt/Mt = μM dt + σM dZt
========================================================================================#

function compute_EψM(x, μx, σx; t::AbstractVector = range(0, 100, step = 1/12), ψ = ones(length(x)), μM = zeros(length(x)), σM = zeros(length(x)))
    feynman_kac_forward(x, μx .+ σM .* σx, σx; t = t, ψ = ψ, V = μM)
end



#========================================================================================

Compute ϵ(x, T) = σD(x) * (σM + σE[M_T | X_t = x])

========================================================================================#

# compute ϵ(x, t) = σD(x) * (σM + σE[M_t | X_0 = x])
function compute_ϵ(x, μx, σx, μM, σM, σD; t::AbstractVector = range(0, 100, step = 1/12))
    u = compute_EψM(x, μx, σx; t = t, μM = μM, σM = σM)
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

