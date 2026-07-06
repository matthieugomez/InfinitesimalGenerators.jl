# Computing expectations (Kolmogorov backward)

Many objects in economics are conditional expectations over the path of a Markov process: forecasts, present values, survival probabilities, option values. The mathematical problem is to compute

```math
u(y, t) = E[\psi(y_T) \mid y_t = y],
```

where ``\psi`` is a terminal payoff or event indicator. This object solves the **Kolmogorov backward equation**

```math
0 = \partial_t u + \mathbb{A} u, \qquad u(\cdot, T) = \psi,
```

where ``\mathbb{A}`` is the infinitesimal generator of the continuous-time process. For instance, if ``y_t`` is the diffusion

```math
dy_t = \mu(y_t) \, dt + \sigma(y_t) \, dZ_t,
```

then the PDE is

```math
0 = \partial_t u + \mu(y) \partial_y u + \frac{1}{2}\sigma(y)^2 \partial_{yy} u,
\qquad u(y, T) = \psi(y).
```

InfinitesimalGenerators.jl solves this PDE on a finite grid. First construct a process object; the package then builds the generator matrix ``\mathbb{A}`` on that grid. For example, take an Ornstein-Uhlenbeck cash flow,

```math
dy_t = \kappa(\bar y - y_t) \, dt + \sigma \, dZ_t,
```

a natural model of a mean-reverting dividend or income stream:

```@example expectations
using InfinitesimalGenerators, LinearAlgebra

κ, σ, ȳ = 0.1, 0.05, 1.0
Y = OrnsteinUhlenbeck(; xbar = ȳ, κ = κ, σ = σ)
ys = only(state_space(Y))
𝔸 = generator(Y)
```

On the grid, ``u_t`` is the vector of values at the grid points, and the backward PDE becomes the linear ODE

```math
0 = \partial_t u_t + \mathbb{A} u_t, \qquad u_T = \psi.
```

Marching backward from ``T`` with implicit Euler gives

```math
u_{t - dt} = (I - dt \, \mathbb{A})^{-1} u_t.
```

The reason to use an implicit step is monotonicity. A finite-difference scheme for a backward equation should preserve the order of payoffs: if ``\psi_1 \leq \psi_2`` at the terminal date, the computed values should satisfy ``u_1 \leq u_2`` at earlier dates. This is the discrete analogue of the maximum principle, and it is one of the standard conditions for convergence of finite-difference schemes for viscosity solutions ([Barles and Souganidis 1991](https://doi.org/10.3233/ASY-1991-4305)).

For a valid Markov generator, ``I - dt \, \mathbb{A}`` is an M-matrix for any ``dt > 0``, so its inverse is non-negative and the implicit scheme is monotone for any time step. By contrast, an explicit step ``u_{t - dt} = (I + dt \, \mathbb{A}) u_t`` is monotone only under the CFL restriction ``dt \leq \min_i 1 / (-\mathbb{A}_{ii})``.

Here is the computation for ``\psi(y) = y``, the conditional mean ``E[y_T \mid y_0]``:

```@example expectations
dt = 0.01
u = collect(ys)
for _ in 1:round(Int, 10.0 / dt)
    u .= (I - dt * 𝔸) \ u
end
nothing # hide
```

For the Ornstein-Uhlenbeck process the conditional mean is known in closed form, ``E[y_T \mid y_0] = \bar y + e^{-\kappa T}(y_0 - \bar y)``, which pins down the accuracy of the discretization:

```@example expectations
maximum(abs, u - (ȳ .+ exp(-κ * 10.0) .* (ys .- ȳ)))
```

The error has two sources: the ``O(dt)`` bias of implicit Euler, which shrinks with the time step, and the reflecting boundaries of the grid, examined below.

## Present values by hand

The same backward equation also handles flow payoffs. The general finite-horizon Feynman-Kac formula is

```math
u(y, t) =
E\left[
\int_t^T e^{-\int_t^s v(y_\tau, \tau) d\tau} f(y_s, s) \, ds
+ e^{-\int_t^T v(y_\tau, \tau) d\tau} \psi(y_T)
\,\Big|\, y_t = y
\right],
```

where ``f`` is a flow payoff, ``v`` is a discount or hazard rate, and ``\psi`` is a terminal payoff. It solves

```math
0 = \partial_t u + \mathbb{A} u - v u + f,
\qquad u(\cdot, T) = \psi.
```

On a finite grid, let ``u_i``, ``f_i``, and ``v_i`` denote the vectors at date ``t_i``, let ``V_i = \operatorname{diag}(v_i)``, and let ``\Delta t_i = t_{i+1} - t_i``. The implicit Euler algorithm is:

1. Set ``u_N = \psi``.
2. For ``i = N-1, \ldots, 0``, solve

```math
u_i = \left(I + \Delta t_i (V_i - \mathbb{A})\right)^{-1} (u_{i+1} + \Delta t_i f_i).
```

The earlier terminal-expectation example is the special case ``f = 0`` and ``v = 0``.

For a time-homogeneous infinite-horizon present value, the algorithm collapses to one resolvent solve:

```math
(V - \mathbb{A}) P = f.
```

With a constant discount rate, ``V = r I``. For example, price a claim to the flow ``y_t`` discounted at rate ``r``:

```math
P(y) = E\left[\int_0^\infty e^{-rt} y_t \, dt \,\Big|\, y_0 = y\right].
```

Differentiating with respect to the starting date gives the stationary backward equation ``r P = y + \mathbb{A} P`` — the continuous-time analogue of "price equals dividend plus discounted expected price". On the grid, solve

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

`feynman_kac` implements the finite-horizon algorithm above. With only a terminal condition, it reproduces the backward march for ``E[y_T \mid y_0]``. The result has one column per date in `ts`, and the first column is the expectation at horizon `ts[end] - ts[1]`:

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

For flow payoffs, pass `f`; for discounting or hazards, pass `v`; for terminal payoffs, pass `ψ`. With `f = ys`, `v = r`, and a horizon long enough that ``e^{-rT} \approx 0``, the helper converges to the resolvent solve above:

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
