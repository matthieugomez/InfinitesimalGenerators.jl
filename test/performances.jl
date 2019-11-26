using InfinitesimalGenerators, Test, Statistics, LinearAlgebra,  Expokit

# all these examples are generators with only real eigenvalues so I should choose better ones

##  Ornstein–Uhlenbeck
κx = 0.1
σ = 0.02
x = range(- 10 * sqrt(σ^2 /(2 * κx)), stop = 10 * sqrt(σ^2 /(2 * κx)), length = 1000)
μx = -κx .* x
σx = σ .* ones(length(x))





## stationnary distribution
@time g = stationary_distribution(x, μx, σx)
#   0.002711 seconds (248 allocations: 713.797 KiB)

## Feynman-Kac
ψ = x.^2
t = range(0, stop = 1000, step = 1/10)
@time u = feynman_kac_forward(x, μx, σx; t = t, ψ = ψ)[:, end]
#   0.019786 seconds (3.09 k allocations: 31.120 MiB, 14.88% gc time)
@test g'u ≈ g'ψ


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
@time l, η, r = principal_eigenvalue(generator_longrun(x, μx, σx, μM, σM; ρ = ρ)(ζ); eigenvector = :both)
#  0.001448 seconds (437 allocations: 660.359 KiB)






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

f = x.^3
𝔸 = generator(x, μx, σx, μM, σM)
t = range(0, stop = 1000, step = 1/10)
u = feynman_kac_forward(𝔸; ψ = f, t = t)[:, end]
l, η, r = principal_eigenvalue(𝔸; eigenvector = :both)
exp(-η * t[end]) * u ./ r
sum(l .* f)
ψhat = stationary_distribution(x, μx .+ σx .* (σM .+ _derive(r, x, μx) ./ r .* σx), σx)
ψtilde = stationary_distribution(x, μx .+ σx .* σM, σx)
l2 = ψhat ./ r ./ sum(ψhat ./ r)
l3 = ψtilde .* r ./ sum(ψtilde .* r)
@test sum(abs2, l - l2) <= 1e-6
@test sum(abs2, l - l3) <= 1e-6 

