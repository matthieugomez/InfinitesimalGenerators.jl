using ContinuousTimeMarkovOperators, Test



# stationnary distribution Ornstein–Uhlenbeck


x = range(-1, 1, length = 1000)
μx = - 0.1 .* x  
σx = 0.2 .* ones(length(x))
stationary_distribution(x, μx, σx)

# Feynman Kac
ψ = x.^2
feynman_kac_forward(x, μx, σx; ψ = ψ)
@time feynman_kac_forward(x, μx, σx; ψ = ψ, t = collect(range(0, 100, step = 1/24)))
feynman_kac_forward(x, μx, σx; ψ = ψ, f = - ones(length(x)))

# One can also compute it using exponential integrator
using Expokit
𝔸 = ContinuousTimeOperators.build_operator(x, zeros(length(x)), μx, 0.5 .* σx.^2)
@time expmv(100.0, 𝔸', ψ)

LawsonEuler(krylov=true, m=50)


feynman_kac_backward(x, μx, σx; ψ = ψ)
feynman_kac_backward(x, μx, σx; ψ = ψ, t = collect(range(0, 100, step = 1/12)))






