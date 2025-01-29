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

#=
Public API
=#

export BoundTuple
export SimulationCondition, planeCouetteFlow, planePoiseuilleFlow

export Transformer
export y2y_F!, y_F2y!
export physical2fourier!, fourier2physical!

export TriDiagonalMatrix, LinearSolver
export solve_system_R!, solve_system_C!

export InputOutputManager
export write_grid
export write_flowfield, read_flowfield!
export write_flow_statistics, write_attribute

export NavierStokesPropagator, TimeStepper
export DomainDescriptor, SimulationCondition, State
export InputOutputManager, Transformer

end
