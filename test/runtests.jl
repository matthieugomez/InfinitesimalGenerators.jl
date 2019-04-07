using InfinitesimalGenerators, Test, BandedMatrices, Statistics, LinearAlgebra,  Expokit


B = BandedMatrix(-1=> 1:5, 2=>1:3)
A = InfinitesimalGenerator(B)
@test isa(A + A, InfinitesimalGenerator)
@test isa(A + I, InfinitesimalGenerator)
@test isa(A * A, BandedMatrix)


##  Ornstein–Uhlenbeck
κx = 0.1
σ = 0.02
x = range(- 3 * sqrt(σ^2 /(2 * κx)), stop = 3 * sqrt(σ^2 /(2 * κx)), length = 1000)
μx = -κx .* x
σx = σ .* ones(length(x))

## stationnary distribution
𝔸 = generator(x, μx, σx)
g = stationary_distribution(𝔸)
@test sum(g .* x) ≈ 0.0 atol = 1e-6
@test sum(g .* x.^2) ≈ σ^2 /(2 * κx) atol = 1e-2

## Feynman-Kac
ψ = x.^2
t = range(0, stop = 100, step = 1/10)
u = feynman_kac_forward(𝔸; t = t, ψ = ψ)
# Check results using exponential integrator. I could also use KrylovKit.exponentiate
𝔸 = generator(x, μx, σx)
@test maximum(abs, u[:, 50] .- expmv(t[50], 𝔸, ψ)) <= 1e-3
@test maximum(abs, u[:, 200] .- expmv(t[200], 𝔸, ψ)) <= 1e-3
@test maximum(abs, u[:, end] .- expmv(t[end], 𝔸, ψ)) <= 1e-5
@test maximum(abs, feynman_kac_forward(𝔸; t = t, ψ = ψ) .- feynman_kac_forward(𝔸; t = collect(t), ψ = ψ)) <= 1e-5


## Multiplicative Functional dM/M = x dt
μM = x
σM = zeros(length(x))
𝔸M = generator(x, μx, σx, μM, σM)
η, f = hansen_scheinkman(𝔸M)
@test η ≈ 0.5 * σ^2 / κx^2 atol = 1e-2
@test maximum(abs, f ./ exp.(x ./ κx) .- mean(f ./ exp.(x ./ κx))) <= 1e-2

t = range(0, stop = 100, step = 1/10)
u = feynman_kac_forward(𝔸M; t = t)
@test log.(stationary_distribution(𝔸)' * u[:, end]) ./ t[end] ≈ η atol = 1e-2

