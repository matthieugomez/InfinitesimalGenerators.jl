# API reference

## Markov processes

```@docs
ContinuousTimeMarkovProcess
ContinuousTimeMarkovChain
DiffusionProcess
OrnsteinUhlenbeck
CoxIngersollRoss
MultivariateDiffusionProcess
ProductProcess
SwitchingProcess
```

## State-space interface

Every process implements the same small shape interface:

- `state_space(X)` returns the state-space axes as a tuple, including for one-dimensional processes.
- `size(X)` returns the tensor shape used for state-shaped arrays.
- `length(X)` returns the number of flattened states, equal to `prod(size(X))`.
- `ndims(X)` returns the number of state-space axes.

```@docs
state_space
```

## The generator and its operators

```@docs
generator
stationary_distribution
feynman_kac
InfinitesimalGenerators.principal_eigenvalue
```

## Additive functionals and tail indices

```@docs
AdditiveFunctionalDiffusion
cgf
tail_index
```

## Finite differences

```@docs
FirstDerivative
SecondDerivative
```
