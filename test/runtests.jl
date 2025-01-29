# Run Test Set

#=
Package
=#

using Test

include("../src/sim/DomainDescriptors.jl")
using ..DomainDescriptors

include("../src/sim/States.jl")
using ..States

include("../src/sim/SimulationConditions.jl")
using ..SimulationConditions

include("../src/util/Transformations.jl")
using ..Transformations

include("../src/util/LinearSolvers.jl")
using ..LinearSolvers

#=
Module Level Unit Test
=#

# DomainDescriptor
include("sim/DomainDescriptors.jl")

# State
include("sim/States.jl")

# Transformer
include("util/Transformations.jl")

# LinearSolver
include("util/LinearSolvers.jl")

#=
Integration Test
=#

include("StateTransformations.jl")
