module InfinitesimalGenerators

using Arpack
using Distributions
using FillArrays
using FiniteDiff
using KrylovKit
using LinearAlgebra
using Roots




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