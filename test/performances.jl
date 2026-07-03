using InfinitesimalGenerators

# all these examples are generators with only real eigenvalues so I should choose better ones

##  Ornstein–Uhlenbeck
κx = 0.1
σ = 0.02
x = range(- 10 * sqrt(σ^2 /(2 * κx)), stop = 10 * sqrt(σ^2 /(2 * κx)), length = 1000)
μx = -κx .* x
σx = σ .* ones(length(x))





## stationnary distribution
@time g = stationary_distribution(DiffusionProcess(x, μx, σx))
#   0.000043 seconds (79 allocations: 193.938 KiB)

## Feynman-Kac
ψ = x.^2
t = range(0, stop = 1000, step = 1/10)
@time u = feynman_kac(generator(DiffusionProcess(x, μx, σx)), t; ψ = ψ)[:, end]
#   0.103 seconds (67 allocations: 76.462 MiB) — scales linearly in length(t);
#   with stop = 100 it runs in 0.0085 seconds (7.8 MiB)
g'u ≈ g'ψ


## test left and right eigenvector
κx = 0.1
σ = 0.02
x = range(- 3 * sqrt(σ^2 /(2 * κx)), stop = 3 * sqrt(σ^2 /(2 * κx)), length = 500)
μx = -κx .* x
σx = σ .* ones(length(x))
μM = -0.01 .+ x
σM = 0.1 .* ones(length(x))
ρ = 1.0
M = AdditiveFunctionalDiffusion(DiffusionProcess(x, μx, σx), μM, σM; ρ = ρ)
@time ζ = tail_index(M)
#  0.001 seconds
@time η, l = cgf_eigenvector(M, ζ, :left)
@time η, r = cgf_eigenvector(M, ζ, :right)
#  0.0001 seconds combined
