# Multivariate processes

Multivariate processes live on tensor-product state spaces. They subtype `ContinuousTimeMarkovProcess{N}`, where `N` is the number of state-space axes: `size(X)` is the shape of the grid (e.g. `(100, 2)` for a diffusion crossed with a two-state chain), and `generator(X)` acts on states flattened in column-major order, so `stationary_distribution(X)` and `feynman_kac(X, ...)` return arrays shaped like `size(X)`.

Start from the univariate building blocks:

```@example multivariate
using InfinitesimalGenerators

X = OrnsteinUhlenbeck(; xbar = 0.0, κ = 0.03, σ = 0.01)
Z = ContinuousTimeMarkovChain([0.5, 1.5], [-0.1 0.1; 0.2 -0.2])
nothing # hide
```

## Products of independent processes

`ProductProcess` combines independent processes into their joint process; the generator is assembled from Kronecker sums, and the state space is the tensor product in argument order:

```@example multivariate
Y = ProductProcess(X, Z)
(ndims(Y), size(Y))
```

## Regime switching

`SwitchingProcess(Z, Xs)` drops the independence: the continuous dynamics *depend* on the chain. In state `only(state_space(Z))[i]`, the process follows `Xs[i]`; all regime processes must share the same grid. This is the natural representation of a solved model with a discrete state — see the [HJB tutorial](hjb.md), where wealth drifts at a policy-implied rate that switches with income:

```@example multivariate
xs = only(state_space(X))
Xlow  = DiffusionProcess(xs, -0.03 .* xs .- 0.01, 0.01 .* ones(length(xs)))
Xhigh = DiffusionProcess(xs, -0.03 .* xs .+ 0.01, 0.01 .* ones(length(xs)))
S = SwitchingProcess(Z, [Xlow, Xhigh])
size(S)
```

## Correlated diffusions

`MultivariateDiffusionProcess` handles correlated diffusions on a tensor-product grid, with drift, variance, and covariance terms already evaluated on it (scalars or arrays), for example after solving an HJB:

```@example multivariate
xs = range(-1, 1, length = 50)
ys = range(0, 2, length = 40)
grid = (; x = xs, y = ys)

μx = repeat(-0.1 .* xs, 1, length(ys))
μy = repeat(reshape(1 .- ys, 1, :), length(xs), 1)
σx = 0.2 .* ones(length(xs), length(ys))
σy = 0.3 .* ones(length(xs), length(ys))
covxy = 0.01 .* ones(length(xs), length(ys))

Y = MultivariateDiffusionProcess(grid;
    drift = (; x = μx, y = μy),
    variance = (; x = σx .^ 2, y = σy .^ 2),
    covariance = (; xy = covxy),
)

g = stationary_distribution(Y)   # an array with size(Y)
size(g)
```

## When is the discretization a valid generator?

`generator(Y)` checks that the discretized operator is a valid Markov generator: rows must sum to zero and off-diagonal entries must be non-negative transition rates. Own-drift and own-variance terms always satisfy this, but the directional cross-derivative stencil used for correlated states is monotone only under a grid-scaled diagonal-dominance condition. For two states `x` and `y`, a useful local rule of thumb is

```julia
abs(covxy) <= variance.x * Δy / Δx
abs(covxy) <= variance.y * Δx / Δy
```

On equally spaced grids this becomes `abs(covxy) <= min(variance.x, variance.y)`, which is stronger than positive semidefiniteness of the covariance matrix. If `generator(Y)` throws a negative-off-diagonal error, the most common fix is to rescale or refine the grids so `Δx / Δy` is closer to `sqrt(variance.x / variance.y)` in the region where the covariance is large. Equivalently, transform states so their local volatilities are more balanced, or reduce the covariance. Use `generator(Y; check = :warn)` or `check = false` only when you intentionally want the raw finite-difference operator without a Markov-process interpretation.

```@example multivariate
@assert abs(sum(g) - 1) <= 1e-8 # hide
nothing # hide
```
