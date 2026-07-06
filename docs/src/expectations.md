# Expected values: the Kolmogorov backward equation

Many objects in economics are conditional expectations over the path of a Markov process: forecasts, present values, survival probabilities, option values. This tutorial computes them *by hand* from the generator matrix — each is a few lines of linear algebra — and then introduces [`feynman_kac`](@ref), the helper function that packages the computation.

Throughout, take a cash flow ``y_t`` that follows an Ornstein–Uhlenbeck process

```math
dy_t = \kappa (\bar y - y_t) \, dt + \sigma \, dZ_t,
```

a natural model of a mean-reverting dividend or income stream. The convenience constructor [`OrnsteinUhlenbeck`](@ref) chooses a grid spanning the stationary distribution and returns the discretized process:

```@example expectations
using InfinitesimalGenerators, LinearAlgebra

κ, σ, ȳ = 0.1, 0.05, 1.0
Y = OrnsteinUhlenbeck(; xbar = ȳ, κ = κ, σ = σ)
ys = only(state_space(Y))
𝔸 = generator(Y)
```

## Forecasts by hand

Fix a horizon ``T`` and consider the forecast ``u(y, t) = E[\psi(y_T) \mid y_t = y]``. It solves the **Kolmogorov backward equation**

```math
0 = \partial_t u + \mathbb{A} u, \qquad u(\cdot, T) = \psi,
```

where, on the grid, ``\mathbb{A}`` is the matrix representation of the infinitesimal operator
``(\mathbb{A}f)(y_i) = \lim_{\Delta t \downarrow 0} \left(E[f(y_{t+\Delta t}) \mid y_t = y_i] - f(y_i)\right) / \Delta t``.
After discretizing the state space, this is a linear ODE in ``\mathbb{R}^n``, solved by marching backward from ``T`` with implicit Euler steps:

```math
u_{t - dt} = (I - dt \, \mathbb{A})^{-1} u_t.
```

In code — factorize once, then one back-substitution per time step:

```@example expectations
function expectation(𝔸, ψ, T; dt = 0.01)
    B = factorize(I - dt * 𝔸)
    u = copy(ψ)
    for _ in 1:round(Int, T / dt)
        u = B \ u
    end
    return u
end

u = expectation(𝔸, collect(ys), 10.0)   # ψ(y) = y: the conditional mean E[y_T | y_0]
nothing # hide
```

For the Ornstein–Uhlenbeck process the conditional mean is known in closed form, ``E[y_T \mid y_0] = \bar y + e^{-\kappa T}(y_0 - \bar y)``, which pins down the accuracy of the discretization:

```@example expectations
maximum(abs, u - (ȳ .+ exp(-κ * 10.0) .* (ys .- ȳ)))
```

The error has two sources: the ``O(dt)`` bias of implicit Euler, which shrinks with the time step, and the reflecting boundaries of the grid, examined below.

## ... and with the helper

`feynman_kac` runs exactly this backward march (reusing the factorization, as above) and returns the whole time path. The terminal condition is the keyword `ψ`; the result has one column per date in `ts`, and the first column is the expectation at horizon `ts[end] - ts[1]`:

```@example expectations
ts = range(0, 10, step = 0.01)
u2 = feynman_kac(Y, ts; ψ = collect(ys))
maximum(abs, u2[:, 1] - u)
```

Because ``\psi`` is arbitrary, probabilities are the same computation — a probability is the expectation of an indicator. The probability that the cash flow is below its long-run mean in ten years:

```@example expectations
prob = feynman_kac(Y, ts; ψ = float.(ys .<= ȳ))[:, 1]
extrema(prob)
```

## Present values by hand

Now price a claim to the flow ``y_t`` discounted at rate ``r``:

```math
P(y) = E\left[\int_0^\infty e^{-rt} y_t \, dt \,\Big|\, y_0 = y\right].
```

Differentiating with respect to the starting date gives the stationary backward equation ``r P = y + \mathbb{A} P`` — the continuous-time analogue of "price equals dividend plus discounted expected price". Discretized, it is a single linear solve in the resolvent of the generator:

```@example expectations
r = 0.05
P = (r * I - 𝔸) \ collect(ys)
nothing # hide
```

The closed form ``P(y) = \bar y / r + (y - \bar y)/(r + \kappa)`` — mean-reverting cash flows are discounted at ``r + \kappa``, not ``r`` — again gives a check, and this time it also reveals a systematic error of the discretization. The package always imposes reflecting boundaries, i.e. a zero derivative at the edges of the grid, while the true ``P`` has slope ``1/(r+\kappa)`` everywhere. The result is a boundary layer: the error is visible at the very edge of the grid, dies out within a few standard deviations, and is negligible where the process actually spends time:

```@example expectations
closedP = ȳ / r .+ (ys .- ȳ) ./ (r + κ)
err = abs.(P - closedP)
interior = abs.(ys .- ȳ) .<= 3 * σ / sqrt(2κ)   # within three sd of the stationary mean
(edge = maximum(err), interior = maximum(err[interior]), average = sum(stationary_distribution(Y) .* err))
```

This is why the grids chosen by `OrnsteinUhlenbeck` and `CoxIngersollRoss` span far into the tails (by default, the ``10^{-10}`` quantiles of the stationary distribution): the boundary distortion then sits where the process essentially never goes.

## ... and with the helper

`feynman_kac` computes the general finite-horizon version,

```math
u(y, 0) = E\left[\int_0^T e^{-\int_0^t v(y_s) ds} f(y_t) \, dt + e^{-\int_0^T v(y_s) ds} \psi(y_T) \,\Big|\, y_0 = y\right],
```

with a flow payoff `f`, a state-dependent (or time-varying) discount rate `v`, and a terminal payoff `ψ`. With `f = ys`, `v = r`, and a horizon long enough that ``e^{-rT} \approx 0``, it converges to the resolvent solve above:

```@example expectations
ts = range(0, 400, step = 0.25)
P2 = feynman_kac(Y, ts; f = collect(ys), v = r .* ones(length(ys)))
maximum(abs, P2[:, 1] - P)
```

The state-dependent discount `v` is what makes the helper more general than the one-line resolvent: it prices claims under a stochastic short rate, computes expected values with state-dependent hazard rates of death or default, and handles time-varying payoffs (pass `f` as a matrix with one column per date).

```@example expectations
@assert maximum(abs, u2[:, 1] - u) <= 1e-10 # hide
@assert maximum(abs, u - (ȳ .+ exp(-κ * 10.0) .* (ys .- ȳ))) <= 1e-2 # hide
@assert maximum(err[interior]) <= 1e-4 # hide
@assert sum(stationary_distribution(Y) .* err) <= 1e-6 # hide
@assert maximum(abs, P2[:, 1] - P) <= 1e-2 # hide
nothing # hide
```
