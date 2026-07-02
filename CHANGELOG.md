# Changelog

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
