# Changelog

## 3.2.0

### Fixed
- `principal_eigenvalue` no longer throws a `SingularException` when the initial shift
  lands exactly on an eigenvalue — which happens whenever the row sums of the Metzler
  matrix are exactly constant, so that the Gershgorin bound used as the initial shift is
  attained (e.g. the tilted generator of a state-independent additive functional over an
  exact chain). The shift is nudged by `tol` and the Rayleigh quotient recovers the exact
  eigenvalue.
- `principal_eigenvalue` on a reducible generator (states that do not communicate) now
  throws an `ArgumentError` explaining that the stationary eigenvector is not unique,
  instead of a bare `SingularException`.
- `feynman_kac` no longer errors when `f` or `v` is a one-column matrix (e.g. a
  state-shaped array with a trailing time dimension of length 1 in the process form);
  a single column is held constant over time.
- `feynman_kac` now throws if `ts` is not increasing; a decreasing grid previously
  flipped implicit Euler into an unstable scheme silently.

### Added
- `tail_index` accepts a `bracket` keyword (default `(1e-5, 1e3)`); if `cgf(m, ξ) - δ`
  has the same sign at both ends, an `ArgumentError` reports the two values instead of a
  cryptic root-finding error.
- The process-level `stationary_distribution(X; ...)` and `feynman_kac(X, ts; ...)`
  forward remaining keyword arguments to `generator(X)`, so `check = :warn` /
  `check = false` reach a `MultivariateDiffusionProcess` without assembling the matrix
  by hand.
- Generator assembly (`generator`, `∂`) now uses the promoted element type of its inputs
  instead of hard-coding `Float64`, so dual numbers flow through: univariate stationary
  distributions, `cgf`, and other spectral objects can be differentiated with
  ForwardDiff with respect to model parameters (multivariate processes assemble
  generically too, but sparse factorization of dual matrices is not supported by the
  ecosystem). Tested against finite differences.
- `cgf` and `cgf_eigenvector` no longer assume that a custom `AdditiveFunctional`
  subtype has an `X` field: the default eigenvector guess is derived from the tilted
  generator, so `tilted_generator` is the only method a custom subtype must define, as
  documented.
- Aqua.jl quality checks (method ambiguities, unbound type parameters, compat bounds,
  stale dependencies) run as part of the test suite; compat entries added for the
  standard libraries and test dependencies.

### Changed
- The zero-row-sum test in `principal_eigenvalue` and the zero-eigenvalue warning in
  `stationary_distribution` now use tolerances scaled by the size of the diagonal (as
  `check_generator` already did), instead of absolute constants — robust to generators
  with very large transition rates (fine grids).
- `feynman_kac` factorizes the implicit-Euler matrix once whenever the time grid is
  uniform and `v` is time-invariant, even when `f` varies over time (previously a
  time-varying `f` forced one factorization per step).
- The positive-semidefiniteness validation in `MultivariateDiffusionProcess` and
  `AdditiveFunctional` reuses a single buffer across grid points and uses a closed form
  instead of `eigmin` for 2×2 problems — faster on large grids, and generic in the
  element type.

## 3.1.0

### Added
- `check_generator(𝔸; atol)` checks that a matrix is a generator (transition-rate)
  matrix — nonnegative off-diagonal entries, rows summing to zero — exploiting matrix
  structure (bands, stored sparse entries) instead of scanning all `n²` entries.
  Violations emit a warning reporting their size rather than throwing, since a matrix
  close to a generator often still yields accurate results. Useful before handing a
  hand-built matrix to `stationary_distribution` or `feynman_kac`;
  `ContinuousTimeMarkovChain` now uses it to check its input (and so warns, rather than
  errors, on an invalid `Q`).
- `SwitchingProcess(Z, Xs)` now accepts *any* process as the modulator `Z`, not just a
  `ContinuousTimeMarkovChain` — e.g. a diffusion modulator yields continuously modulated
  dynamics, equivalent to a `MultivariateDiffusionProcess` with independent innovations.
- `AdditiveFunctional` now works for *any* `ContinuousTimeMarkovProcess` — chains,
  products, switching processes, multivariate diffusions — not just univariate
  diffusions, and gains a keyword constructor with canonical coefficients:
  `AdditiveFunctional(X; drift, variance, covariance)`, where `drift` and `variance`
  are scalars or arrays shaped like `size(X)`. The `covariance` keyword specifies
  `cov(dm, dx)` directly and — new capability — also works for
  `MultivariateDiffusionProcess` states (a `NamedTuple` such as `(; x = cx)`).
  The positional `AdditiveFunctional(X, μm, σm; ρ)` remains as SDE-style sugar.
- `tilted_generator(m, ξ)` is exported: the tilted generator matrix `𝔸_ξ` from which
  all additive-functional operators are computed, and the extension point for custom
  functionals (parallel to `generator` for custom processes). It can be passed
  directly to `feynman_kac` for finite-horizon moments `E[e^{ξ mₜ} ψ(xₜ)]`.
- `principal_eigenvalue` is exported (previously documented under its qualified
  name). It now checks on entry that its input is a Metzler matrix (nonnegative
  off-diagonal entries) — the Perron–Frobenius hypothesis its result relies on —
  warning if not.
- The `direction` keyword of `FirstDerivative` and `SecondDerivative` accepts
  `:forward`/`:up` and `:backward`/`:down` as synonyms everywhere (the `:up`/`:down`
  vocabulary matches EconPDEs' `_up`/`_down` derivative names).

### Changed
- `cgf(m, ξ)` now returns the scalar `Λ(ξ)`. The old closure form `cgf(m)(ξ)`,
  returning a tuple, is deprecated; eigenvectors are available through the new
  `cgf_eigenvector(m, ξ, :right)` / `cgf_eigenvector(m, ξ, :left)` (the left
  eigenvector is normalized to sum to one, as before).
- `AdditiveFunctionalDiffusion` is a legacy type: `AdditiveFunctional` constructs an
  equivalent functional for diffusion states, including correlated noise. The type
  remains exported and functional.
- The rebirth-distribution keyword of `stationary_distribution` is now `rebirth`
  instead of `ψ`, to avoid the collision with `feynman_kac`, where `ψ` is the
  terminal payoff. The old keyword still works with a deprecation warning.
- The multi-dimensional methods of `FirstDerivative` and `SecondDerivative` now
  return lazy arrays (entries computed on demand), like the one-dimensional
  methods, instead of materialized `Array`s. Use `collect` to materialize.
- Cross-derivatives (`SecondDerivative(grid, F, dim1, dim2)` with `dim1 != dim2`)
  now throw if a nonzero `bc` is passed; the keyword was previously accepted and
  silently ignored.

## 3.0.0

### Breaking
- `state_space(X)` now always returns a tuple of state-space axes, including for
  one-dimensional processes. Use `only(state_space(X))` to recover the grid of a
  one-dimensional process. `size(X)` is now derived once from `state_space(X)`.
- `ContinuousTimeMarkovChain` now stores its state labels in `states` instead
  of `z`.
- `DiffusionProcess` and `AdditiveFunctionalDiffusion` are now immutable
  `struct`s (previously `mutable struct`), so the invariants checked in their
  constructors cannot be bypassed. Code that reassigned their fields after
  construction will no longer work.
- `ProductProcess` now represents the independent product of Markov processes
  in argument order, e.g. `ProductProcess(X, Z)`.
- `generator(X::MultivariateDiffusionProcess)` now throws by default when the
  correlated-state stencil creates negative off-diagonal rates. Pass
  `check = :warn` or `check = false` to inspect the raw finite-difference
  operator anyway.
- The abstract process hierarchy is now `ContinuousTimeMarkovProcess{N}`,
  where `N` is the number of tensor-product state-space axes. The old
  `MarkovProcess`, `UnivariateMarkovProcess`, and `MultivariateMarkovProcess`
  abstract names were removed.

### Added
- `ContinuousTimeMarkovProcess{N}` as the dimension-parametric abstract
  supertype for process dispatch.
- `ContinuousTimeMarkovChain` as the explicit finite-state continuous-time
  chain type.
- `MultivariateDiffusionProcess` for tensor-product diffusion grids with drift,
  variance, and covariance arrays.
- `FirstDerivative((xs, ys), F, 1)` and
  `SecondDerivative((xs, ys), F, 1, 2)` forms for multidimensional finite
  differences on tuple grids. `NamedTuple` grids with symbolic dimensions remain
  supported as convenience syntax.
- Process-level `feynman_kac(X::ContinuousTimeMarkovProcess, ts; ...)`,
  accepting state-shaped arrays and returning state-shaped time paths.

### Fixed
- `∂(::DiffusionProcess)` no longer produces `NaN` rows at nodes where the drift
  is exactly zero. The operator is now built directly (rather than as
  `Diagonal(μx) \ generator(…)`, which divided by zero there): it uses an upwind
  scheme matching `generator`, and falls back to a central difference at interior
  zero-drift nodes.
- One-dimensional diffusion generators now consistently drop outward drift at
  reflecting boundaries, matching the multidimensional generator.
- `MultivariateDiffusionProcess` construction now always goes through validation
  and rejects non-PSD covariance matrices.
- CI: the nightly Julia job is now correctly allowed to fail
  (`continue-on-error`) instead of failing the whole workflow. Previously the
  `allow_failure` matrix key was defined but never referenced.

### Changed
- `feynman_kac` allocates its output using the promoted element type of its
  inputs, so `Float32` (and other) element types are preserved instead of being
  forced to `Float64`.
- `stationary_distribution(X::ContinuousTimeMarkovProcess)` now returns arrays
  with `size(X)` for multidimensional processes, while matrix-level calls remain
  flat.
- The test suite is organized into `@testset`s and tightened: several
  expressions that looked like assertions but tested nothing are now real
  `@test`s, and the zero-drift `∂` path and `feynman_kac` element type are
  covered.
