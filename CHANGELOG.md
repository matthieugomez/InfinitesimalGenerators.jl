# Changelog

## 2.5.2

### Fixed
- `∂(::DiffusionProcess)` no longer produces `NaN` rows at nodes where the drift
  is exactly zero. The operator is now built directly (rather than as
  `Diagonal(μx) \ generator(…)`, which divided by zero there): it uses an upwind
  scheme matching `generator`, and falls back to a central difference at interior
  zero-drift nodes.
- CI: the nightly Julia job is now correctly allowed to fail
  (`continue-on-error`) instead of failing the whole workflow. Previously the
  `allow_failure` matrix key was defined but never referenced.

### Changed
- `feynman_kac` allocates its output using the promoted element type of its
  inputs, so `Float32` (and other) element types are preserved instead of being
  forced to `Float64`.
- `DiffusionProcess` and `AdditiveFunctionalDiffusion` are now immutable
  `struct`s, so the invariants checked in their constructors cannot be bypassed.
- The test suite is organized into `@testset`s and tightened: several
  expressions that looked like assertions but tested nothing are now real
  `@test`s, and the zero-drift `∂` path and `feynman_kac` element type are
  covered.
