# Run Test Set

#=
Package
=#

using Test

include("../src/sim/DomainDescriptors.jl")
using ..DomainDescriptors

include("../src/sim/States.jl")
using ..States

include("../src/util/Transformations.jl")
using ..Transformations

#=
Module Level Unit Test
=#

# DomainDescriptor
include("sim/DomainDescriptors.jl")

# State
include("sim/States.jl")

# Transformer
include("util/Transformations.jl")

#=
Integration Test
=#

include("StateTransformations.jl")
