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

Compute Hansen Scheinkmann decomposition M = e^{ηt}f(x_t)W_t
Return f, η, g

========================================================================================#
# Compute η in Hansen Scheinkmann decomposition M = e^{ηt}f(x_t)W_t
function compute_η(x, μx, σx, μM, σM; method = :krylov, eigenvector = :right)
    n = length(x)
    T = zeros(n, n)
    Δ = EconPDEs.make_Δ(x)
    compute_η!(T, Δ, μx, σx, μM, σM; method = method, eigenvector = eigenvector)
end

function compute_η!(T, Δ, μx, σx, μM, σM; method = :krylov, eigenvector = :right)
    build_operator!(T, Δ, μM, σM .* σx .+ μx, 0.5 .* σx.^2)
    principal_eigenvalue(T; method = method, eigenvector = eigenvector)
end

#========================================================================================

Compute ϵ(x, T) = σD(x) * (σM + σE[M_T | X_t = x])

========================================================================================#


# compute ϵ(x, t) = σD(x) * (σM + σE[M_t | X_0 = x])
function compute_ϵ(x, μx, σx, μM, σM, σD; Y = 100, P = 4)
    t, E = compute_EψM(x, μx, σx; μM = μM, σM = σM, Y = Y, P = P)
    for i in 1:(Y * P)
        E[:, i] = σD .* (σM .+ _derive(E[:, i], x, μx) ./ E[:, i] .* σx)
    end
    return t, E
end

