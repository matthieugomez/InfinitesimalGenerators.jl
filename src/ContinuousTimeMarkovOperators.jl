module ContinuousTimeMarkovOperators
using LinearAlgebra, BandedMatrices, KrylovKit

##############################################################################
##
## Load files
##
##############################################################################
include("build_operator.jl")
include("stationary_distribution.jl")
include("feynman_kac.jl")
include("geometric_functionals.jl")


##############################################################################
##
## Exported methods and types 
##
##############################################################################
export stationary_distribution,
feynman_kac_backward,
feynman_kac_forward,
compute_EψM,
compute_η,
compute_ϵ 
end