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
Module Level Test
=#

# DomainDescriptor
include("sim/TestDomainDescriptors.jl")

# State
include("sim/TestStates.jl")

# Transformer
include("util/TestTransformations.jl")

#=
Integration Test
=#

