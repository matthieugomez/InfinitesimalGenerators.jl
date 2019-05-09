using InfinitesimalGenerators, Test, Statistics, LinearAlgebra,  Expokit


##  Ornstein–Uhlenbeck
κx = 0.1
σ = 0.02
x = range(- 10 * sqrt(σ^2 /(2 * κx)), stop = 10 * sqrt(σ^2 /(2 * κx)), length = 1_000)
μx = -κx .* x
σx = σ .* ones(length(x))





## stationnary distribution
@time g = stationary_distribution(x, μx, σx)
#   0.002711 seconds (248 allocations: 713.797 KiB)

## Feynman-Kac
ψ = x.^2
t = range(0, stop = 100, step = 1/10)
@time u = feynman_kac_forward(x, μx, σx; t = t, ψ = ψ)
#   0.019786 seconds (3.09 k allocations: 31.120 MiB, 14.88% gc time)


## test left and right eigenvector
κx = 0.1
σ = 0.02
x = range(- 3 * sqrt(σ^2 /(2 * κx)), stop = 3 * sqrt(σ^2 /(2 * κx)), length = 500)
μx = -κx .* x
σx = σ .* ones(length(x))
μM = -0.01 .+ x
σM = 0.1 .* ones(length(x))
ρ = 1.0
@time ζ = tail_index(x, μx, σx, μM, σM; ρ = ρ)
#  0.094506 seconds (15.21 k allocations: 24.270 MiB, 40.79% gc time)
@time g, η, f = principal_eigenvalue(generator_longrun(x, μx, σx, μM, σM; ρ = ρ)(ζ); eigenvector = :both)
#  0.001448 seconds (437 allocations: 660.359 KiB)




