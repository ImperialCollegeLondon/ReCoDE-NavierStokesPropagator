module NavierStokesPropagator

#=
Package
=#

include("sim/DomainDescriptors.jl")
using .DomainDescriptors

include("sim/States.jl")
using .States

include("sim/SimulationConditions.jl")
using .SimulationConditions

# include("util/InputOutputManagers.jl")
# using .InputOutputManagers

include("util/Transformations.jl")
using .Transformations

include("util/LinearSolvers.jl")
using .LinearSolvers

# include("util/IntegrationSuite.jl")
# using .IntegrationSuite

#=
Public API
=#

export DomainDescriptor

export State

export SimulationCondition
export BoundTuple
export SimulationCondition, planeCouetteFlow, planePoiseuilleFlow
export init_flowfield!, apply_perturbation!

export Transformer
export y2y_F!, y_F2y!
export physical2fourier!, fourier2physical!

export LinearSolvers
export TriDiagonalMatrix, LinearSolver
export define_sys_inter_v!, solve_sys_inter_v!
export define_sys_inter_uw!, solve_sys_inter_uw!
export define_sys_pressure!, solve_sys_pressure_update!, solve_sys_pressure_poisson!

export NavierStokesPropagator

end
