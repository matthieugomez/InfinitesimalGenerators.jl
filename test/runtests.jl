using InfinitesimalGenerators, Test, Statistics, LinearAlgebra,  Expokit


κ = 0.1
σ = 0.02
X = OrnsteinUhlenbeck(;κ =κ, σ = σ, length = 1000)


## Feynman-Kac
ψ = X.x.^2
t = range(0, stop = 100, step = 1/10)
u = feynman_kac(X; t = t, ψ = ψ, direction = :forward)
@test maximum(abs, u[:, 50] .- expmv(t[50], generator(X), ψ)) <= 1e-3
@test maximum(abs, u[:, 200] .- expmv(t[200], generator(X), ψ)) <= 1e-3
@test maximum(abs, u[:, end] .- expmv(t[end], generator(X), ψ)) <= 1e-5
@test maximum(abs, feynman_kac(X; t = t, ψ = ψ, direction = :forward) .- feynman_kac(X; t = collect(t), ψ = ψ, direction = :forward)) <= 1e-5


## Multiplicative Functional dM/M = x dt
M = MultiplicativeFunctional(X, X.x, zeros(length(X.x)))
l, η, r = hansen_scheinkman(M)
@test η ≈ 0.5 * σ^2 / κ^2 atol = 1e-2
r_analytic = exp.(X.x ./ κ) 
@test norm(r ./ sum(r) .- r_analytic ./ sum(r_analytic)) <= 2 * 1e-3
t = range(0, stop = 200, step = 1/10)
u = feynman_kac(M; t = t, direction = :forward)
@test log.(stationary_distribution(X)' * u[:, end]) ./ t[end] ≈ η atol = 1e-2


## test left and right eigenvector with correlation
M = MultiplicativeFunctional(X, X.x, 0.01 * ones(length(X.x)); ρ = 1)
l, η, r = hansen_scheinkman(M; eigenvector = :both)
ψ_tilde = stationary_distribution(MarkovProcess(X.x, X.μx .+ M.ρ .* M.σM .* X.σx, X.σx))
@test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3


# Test Multiplicative with ρ = 0.0
μM = -0.01
σM = 0.1
M = MultiplicativeFunctional(X, μM .+ X.x, σM .* ones(length(X.x)))
ζ = tail_index(M)
ζ_analytic = 2 * (-μM + σM^2/2) / (σM^2 + (σ / κ)^2)
@test ζ ≈ ζ_analytic atol = 1e-2
l, η, r = cgf_longrun(M; eigenvector = :both)(ζ)
@test η ≈ 0.0 atol = 1e-5
ψ = stationary_distribution(X)
@test (r .* ψ) ./ sum(r .* ψ) ≈ l rtol = 1e-3



# Test that the modified process μ + σ^2 ∂ ln(r) has a stationary distribution given by $r^2ψ$
X = OrnsteinUhlenbeck(;κ =κ, σ = σ, length = 1000)
M = MultiplicativeFunctional(X, μM .+ X.x .- 0.02, σM .* ones(length(X.x)))
ψ = stationary_distribution(X)
ζ = tail_index(M)
l, η, r = cgf_longrun(M; eigenvector = :both)(ζ)
ψ_cond = stationary_distribution(MarkovProcess(X.x, X.μx .+ X.σx.^2 .* InfinitesimalGenerators.∂(X, log.(r)), X.σx))
@test (r.^2 .* ψ) ./ sum(r.^2 .* ψ) ≈ ψ_cond rtol = 1e-1


# Test Multiplicative with ρ = 1.0
M = MultiplicativeFunctional(X, M.μM, M.σM; ρ = 1.0)
ζ = tail_index(M)
l, η, r = cgf_longrun(M; eigenvector = :both)(ζ)
@test η ≈ 0.0 atol = 1e-5
ψ_tilde = stationary_distribution(MarkovProcess(X.x, X.μx .+ ζ .* M.σM .* M.ρ .* X.σx , X.σx))
@test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3

# Test CIR
gbar = 0.03
σ = 0.01
X = CoxIngersollRoss(xbar = gbar, κ = κ, σ = σ)
M = MultiplicativeFunctional(X, X.x, zeros(length(X)))
η_analytic = gbar * κ^2 / σ^2 * (1 - sqrt(1 - 2 * σ^2 / κ^2))
@test cgf_longrun(M)(1.0)[2] ≈ η_analytic rtol = 1e-2


# for CIR the speed is given by 
# speed_analytic = - (g .- 0.009) + xbar * κ / sqrt(κ^2 - 2 * σ^2 * ζ)
# r, η, l = cgf_longrun(M, eigenvector = :both)(ζ)