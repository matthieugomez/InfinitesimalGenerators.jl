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

```@docs
state_space
Base.size(::ContinuousTimeMarkovProcess)
Base.length(::ContinuousTimeMarkovProcess)
Base.ndims(::ContinuousTimeMarkovProcess)
```

## The generator and its operators

```@docs
generator
check_generator
stationary_distribution
feynman_kac
jointoperator
principal_eigenvalue
```

## Additive functionals and tail indices

```@docs
AdditiveFunctional
tilted_generator
cgf
cgf_eigenvector
tail_index
AdditiveFunctionalDiffusion
```

## Finite differences

```@docs
FirstDerivative
SecondDerivative
```
