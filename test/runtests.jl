using ContinuousTimeOperators, Test



# stationnary distribution Ornstein–Uhlenbeck
x = range(-1, 1, length = 100)
μx = - 0.1 .* x  
σx = 0.2 .* ones(length(x))
stationary_distribution(x, μx, σx)


# Feynman Kac
ψ = x.^2
feynman_kac_backward(x, μx, σx; ψ = ψ)
feynman_kac_forward(x, μx, σx; ψ = ψ)



