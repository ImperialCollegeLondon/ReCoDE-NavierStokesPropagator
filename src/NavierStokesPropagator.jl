module NavierStokesPropagator

#=
Package
=#

include("sim/DomainDescriptors.jl")
using .DomainDescriptors

include("sim/States.jl")
using .States

# include("sim/SimulationConditions.jl")
# using .SimulationConditions

# include("util/InputOutputManagers.jl")
# using .InputOutputManagers

include("util/Transformations.jl")
using .Transformations

# include("util/LinearSolvers.jl")
# using .LinearSolvers

# include("util/IntegrationSuite.jl")
# using .IntegrationSuite

using Random

#=
Public API
=#

export DomainDescriptor

export State

export Transformer
export y2y_F!, y_F2y!
export physical2fourier!, fourier2physical!

export NavierStokesPropagator

end
