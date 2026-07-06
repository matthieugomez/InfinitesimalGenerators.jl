# Distribution dynamics: the Kolmogorov forward equation

The previous tutorial computed expectations, which march *backward* in time. Cross-sectional distributions march *forward*: given today's distribution of income, wealth, or firm size across agents, where is the distribution next year, and where does it settle in the long run? This tutorial computes both by hand from the generator matrix, and then introduces the helper [`stationary_distribution`](@ref).

Take the same Ornstein–Uhlenbeck process as before — now interpreted as the log productivity of a cross-section of firms:

```@example distributions
using InfinitesimalGenerators, LinearAlgebra

κ, σ = 0.1, 0.05
X = OrnsteinUhlenbeck(; xbar = 0.0, κ = κ, σ = σ)
xs = only(state_space(X))
𝔸 = generator(X)
```

## Evolving a distribution by hand

If expectations satisfy the backward equation ``\partial_t u = \mathbb{A} u``, densities satisfy its adjoint, the **Kolmogorov forward (Fokker–Planck) equation**:

```math
\partial_t g = \mathbb{A}' g.
```

The transpose is not incidental: ``E[f(x_t)] = g' f`` for any test function, so whatever operator advances expectations backward, its transpose advances the distribution forward. On the discretized state space, `g` is the vector of probability masses at each grid point, and the forward equation is again a linear ODE solved with implicit Euler steps:

```@example distributions
function evolve(𝔸, g0, T; dt = 0.1)
    B = factorize(I - dt * copy(𝔸'))
    g = copy(g0)
    for _ in 1:round(Int, T / dt)
        g = B \ g
    end
    return g
end
nothing # hide
```

Start every firm at the same productivity — a point mass one standard deviation above the mean — and watch the cross-section spread out and recenter:

```@example distributions
using Plots

g0 = zeros(length(xs))
g0[findmin(abs.(xs .- σ / sqrt(2κ)))[2]] = 1.0

Δx = step(xs)
plot(xlabel = "log productivity x", ylabel = "cross-sectional density")
for T in (1, 5, 20, 100)
    plot!(xs, evolve(𝔸, g0, T) ./ Δx; label = "t = $T")
end
current()
```

The division by `Δx` converts masses to a *density* for plotting: `g` sums to one without grid weights, so the density at a grid point is the mass divided by the cell width. On this uniform grid the two differ only by a constant factor; on a non-uniform grid (as in the [HJB tutorial](hjb.md)) the distinction matters, since raw masses would trace the grid spacing rather than the shape of the distribution.

The discretized process is an honest Markov chain — rows of ``\mathbb{A}`` sum to zero and off-diagonals are non-negative — so the masses stay non-negative and sum to one at every date, with no renormalization needed:

```@example distributions
sum(evolve(𝔸, g0, 100))
```

For the Ornstein–Uhlenbeck process the transition distribution is known in closed form — Gaussian with mean ``e^{-\kappa t} x_0`` and variance ``\frac{\sigma^2}{2\kappa}(1 - e^{-2\kappa t})`` — which checks the discretization at, say, ``t = 5``:

```@example distributions
t, x0 = 5.0, xs[findmax(g0)[2]]
gt = evolve(𝔸, g0, t)
mean_t = sum(gt .* xs)
var_t = sum(gt .* xs .^ 2) - mean_t^2
(; mean_t, closed_mean = exp(-κ * t) * x0,
   var_t, closed_var = σ^2 / (2κ) * (1 - exp(-2κ * t)))
```

The mean is essentially exact: the upwind scheme discretizes the drift ``\kappa(\bar x - x)`` without bias. The variance is a few percent too high — upwinding adds a little numerical diffusion, the price paid for a discretization that is guaranteed to be a well-defined Markov chain (masses stay non-negative no matter how coarse the grid). Refining the grid shrinks the gap.

## The stationary distribution by hand

As ``t \to \infty`` the distribution converges to the stationary distribution, the fixed point of the forward equation:

```math
\mathbb{A}' g = 0, \qquad \textstyle\sum_i g_i = 1.
```

That is a linear system: `g` is the null vector of the transposed generator, normalized to sum to one. By hand, replace one (redundant) row of ``\mathbb{A}'`` with the normalization:

```@example distributions
B = Matrix(𝔸')
B[end, :] .= 1.0
g∞ = B \ [zeros(length(xs) - 1); 1.0]
nothing # hide
```

For the Ornstein–Uhlenbeck process this should be the ``N(0, \sigma^2/2\kappa)`` density.
The maximum gap, relative to the density's peak:

```@example distributions
s2∞ = σ^2 / (2κ)
closed∞ = @. exp(-xs^2 / (2s2∞)) / sqrt(2π * s2∞)
maximum(abs, g∞ ./ Δx - closed∞) / maximum(closed∞)
```

The same couple of percent as the variance check above — upwinding's numerical diffusion —
and it halves with each doubling of the grid.

## ... and with the helper

`stationary_distribution` performs this computation for any process in the package (using a robust eigenvalue method rather than the row-replacement trick, and reshaping the result to the shape of the state space for multivariate processes):

```@example distributions
g = stationary_distribution(X)
maximum(abs, g - g∞)
```

Two variations are worth knowing:

- **Convergence check.** The evolved distribution approaches the stationary one at the speed of mean reversion: `maximum(abs, evolve(𝔸, g0, 100) - g)` is of the order of the autocorrelation ``e^{-\kappa \cdot 100} \approx 5 \times 10^{-5}``, and by ``t = 300`` the two agree to machine precision.
- **Death and rebirth.** `stationary_distribution(X; δ = δ, rebirth = ψ)` computes the stationary distribution when agents die at rate ``\delta`` and are reborn with distribution ``\psi`` — the resolvent ``(\delta I - \mathbb{A}')^{-1} \delta \psi`` — which is the relevant object in perpetual-youth and firm-entry models.

```@example distributions
@assert abs(sum(evolve(𝔸, g0, 100)) - 1) <= 1e-10 # hide
@assert maximum(abs, g - g∞) <= 1e-10 # hide
@assert abs(mean_t - exp(-κ * t) * x0) <= 1e-3 # hide
@assert abs(var_t / (σ^2 / (2κ) * (1 - exp(-2κ * t))) - 1) <= 0.1 # hide
@assert maximum(abs, evolve(𝔸, g0, 100) - g) <= 1e-5 # hide
@assert maximum(abs, evolve(𝔸, g0, 300) - g) <= 1e-12 # hide
nothing # hide
```
