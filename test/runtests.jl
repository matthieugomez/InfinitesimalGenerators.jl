using LinearAlgebra, ContinuousTimeMarkovOperators, Expokit, Test


##  Ornstein–Uhlenbeck
κx = 0.1
σ = 0.2
x = range(- 3 * sqrt(σ^2 /(2 * κx)), stop = 3 * sqrt(σ^2 /(2 * κx)), length = 1000)
μx = -κx .* x
σx = σ .* ones(length(x))

## stationnary distribution
g = stationary_distribution(x, μx, σx)
@test sum(g .* x) ≈ 0.0 atol = 1e-6
@test sum(g .* x.^2) ≈ σ^2 /(2 * κx) atol = 1e-2

## Feynman-Kac
ψ = x.^2
t = range(0, 100, step = 1/100)
u = feynman_kac_forward(x, μx, σx; ψ = ψ, t = t)
# Check results using exponential integrator
𝔸 = generator(x, μx, σx)
@test maximum(abs, u[:, 50] .- expmv(t[50], 𝔸, ψ)) <= 1e-3
@test maximum(abs, u[:, 200] .- expmv(t[200], 𝔸, ψ)) <= 1e-3
@test maximum(abs, u[:, end] .- expmv(t[end], 𝔸, ψ)) <= 1e-5
@test maximum(abs, feynman_kac_forward(x, μx, σx; ψ = ψ, t = t) .- feynman_kac_forward(x, μx, σx; ψ = ψ, t = collect(t))) <= 1e-5








