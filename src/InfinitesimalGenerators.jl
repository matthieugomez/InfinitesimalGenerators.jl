module InfinitesimalGenerators

using Arpack: eigs
using Distributions: Normal, Gamma, quantile
using FillArrays: Ones, Zeros
using KrylovKit: eigsolve
using LinearAlgebra: Diagonal, Tridiagonal, I, factorize, ldiv!, diag
using Roots: fzero




include("MarkovProcess.jl")
include("AdditiveFunctional.jl")
include("feynman_kac.jl")
include("principal_eigenvalue.jl")



export 
MarkovProcess,
generator,
state_space,
stationary_distribution,
feynman_kac,
DiffusionProcess,
OrnsteinUhlenbeck,
CoxIngersollRoss,
AdditiveFunctional,
cgf,
tail_index,
AdditiveFunctionalDiffusion
end