using InfinitesimalGenerators, Test, Statistics, LinearAlgebra,  Expokit

xbar = 0.0
κ = 0.1
σ = 0.02
X = OrnsteinUhlenbeck(; xbar = xbar, κ = κ, σ = σ, length = 1000)


## Feynman-Kac
ψ = X.x.^2
t = range(0, stop = 100, step = 1/10)
u = feynman_kac(generator(X); t = t, ψ = ψ, direction = :forward)
@test maximum(abs, u[:, 50] .- expmv(t[50], generator(X), ψ)) <= 1e-3
@test maximum(abs, u[:, 200] .- expmv(t[200], generator(X), ψ)) <= 1e-3
@test maximum(abs, u[:, end] .- expmv(t[end], generator(X), ψ)) <= 1e-5
@test maximum(abs, feynman_kac(X; t = t, ψ = ψ, direction = :forward) .- feynman_kac(X; t = collect(t), ψ = ψ, direction = :forward)) <= 1e-5


## Multiplicative Functional dM/M = x dt
m = AdditiveFunctionalDiffusion(X, X.x, zeros(length(X.x)))
l, η, r = cgf(m; eigenvector = :right)(1)
@test η ≈ xbar + 0.5 * σ^2 / κ^2 atol = 1e-2
r_analytic = exp.(X.x ./ κ) 
@test norm(r ./ sum(r) .- r_analytic ./ sum(r_analytic)) <= 2 * 1e-3
t = range(0, stop = 200, step = 1/10)
u = feynman_kac(m; t = t, direction = :forward)(1)
@test log.(stationary_distribution(X)' * u[:, end]) ./ t[end] ≈ η atol = 1e-2


# test speed
μm = - 0.06
δ = 0.0
m = AdditiveFunctionalDiffusion(X, X.x .+ μm, zeros(length(X.x)), δ = δ)
ζ = tail_index(m)
@test μm * ζ + 0.5 * ζ^2 * (σ^2 / κ^2) - δ ≈ 0.0 atol = 1e-2
l, _, r = cgf(m; eigenvector = :both)(ζ)
f =  exp.(ζ .* X.x ./ κ)
norm(f ./ sum(f) .- r ./ sum(r)) <= 1e-2
ψ_reaching = r.* l ./ sum(r .* l)
speed = sum(ψ_reaching .* m.μm) 
@test speed  ≈ μm  +  ζ * (σ^2 / κ^2) atol = 1e-2


# test transformation with a funciton # m2 = p * M. 
# This does not work very well. Note that it works only if the distribution of p has a thinner tail than the distirbuiton of M
m = AdditiveFunctionalDiffusion(X, X.x .- 0.06, zeros(length(X.x)))
ζ = tail_index(m)
l, η, r = cgf(m; eigenvector = :both)(ζ)
p =  exp.(5 .* X.x)
m2 = AdditiveFunctionalDiffusion(X, m.μm .+ (generator(m.X) * log.(p)),  (InfinitesimalGenerators.∂(X) * log.(p)) .* X.σx, ρ = 1)
l2, η2, r2 = cgf(m2; eigenvector = :both)(ζ)
r3 = (r ./ p.^ζ) ./ sum(r ./ p.^ζ)
r2 = r2 ./ sum(r2)
#@test r2 ≈ r3 rtol = 1e-2
l3 = (l .* p.^ζ) ./ sum(l .* p.^ζ)
#@test l2 ≈ l3 rtol = 1e-1

l' * (generator(m.X) * log.(p))
ψ_reaching = l.*r ./ sum(l.* r)
ψ_reaching' * (generator(m.X) * log.(p))

## test left and right eigenvector with correlation
m = AdditiveFunctionalDiffusion(X, X.x, 0.01 * ones(length(X.x)); ρ = 1)
l, η, r = cgf(m; eigenvector = :both)(1)
ψ_tilde = stationary_distribution(DiffusionProcess(X.x, X.μx .+ m.ρ .* m.σm .* X.σx, X.σx))
@test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3


# Test Multiplicative with ρ = 0.0
μm = -0.01
σm = 0.1
m = AdditiveFunctionalDiffusion(X, μm .+ X.x, σm .* ones(length(X.x)))
ζ = tail_index(m)
ζ_analytic = 2 * (-μm) / (σm^2 + (σ / κ)^2)
@test ζ ≈ ζ_analytic atol = 1e-2
l, η, r = cgf(m; eigenvector = :both)(ζ)
@test η ≈ 0.0 atol = 1e-4
ψ = stationary_distribution(X)
@test (r .* ψ) ./ sum(r .* ψ) ≈ l rtol = 1e-3



# Test that the modified process μ + σ^2 ∂ ln(r) has a stationary distribution given by $r^2ψ$
X = OrnsteinUhlenbeck(;κ =κ, σ = σ, length = 1000)
m = AdditiveFunctionalDiffusion(X, μm .+ X.x .- 0.02, σm .* ones(length(X.x)))
ψ = stationary_distribution(X)
ζ = tail_index(m)
l, η, r = cgf(m; eigenvector = :both)(ζ)
ψ_cond = stationary_distribution(DiffusionProcess(X.x, X.μx .+ X.σx.^2 .* (InfinitesimalGenerators.∂(X) * log.(r)), X.σx))
@test (r.^2 .* ψ) ./ sum(r.^2 .* ψ) ≈ ψ_cond rtol = 1e-1


# Test Multiplicative with ρ = 1.0
m = AdditiveFunctionalDiffusion(X, m.μm, m.σm; ρ = 1.0)
ζ = tail_index(m)
l, η, r = cgf(m; eigenvector = :both)(ζ)
@test η ≈ 0.0 atol = 1e-3
ψ_tilde = stationary_distribution(DiffusionProcess(X.x, X.μx .+ ζ .* m.σm .* m.ρ .* X.σx , X.σx))
@test (r .* ψ_tilde) ./ sum(r .* ψ_tilde) ≈ l rtol = 1e-3

# Test CIR
gbar = 0.03
σ = 0.01
X = CoxIngersollRoss(xbar = gbar, κ = κ, σ = σ)
m = AdditiveFunctionalDiffusion(X, X.x, zeros(length(X.x)))
η_analytic = gbar * κ^2 / σ^2 * (1 - sqrt(1 - 2 * σ^2 / κ^2))
@test cgf(m)(1.0)[2] ≈ η_analytic rtol = 1e-2


# for CIR the speed is given by 
# speed_analytic = - (g .- 0.009) + xbar * κ / sqrt(κ^2 - 2 * σ^2 * ζ)
# r, η, l = cgf(m, eigenvector = :both)(ζ)