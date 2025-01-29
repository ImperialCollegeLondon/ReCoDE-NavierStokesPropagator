# Documenter make.jl

#=
Package
=#

include("../src/NavierStokesPropagator.jl")
using .NavierStokesPropagator

using Documenter

#=
Make Documenter
=#

makedocs(sitename="Navier-Stokes Propagator", pages=[
    "Home" => "index.md",
    "Simulation Constructs" => Any["sim/DomainDescriptors.md", "sim/States.md", "sim/SimulationConditions.md"],
    "Simulation Utilities" => Any["util/Transformations.md", "util/LinearSolvers.md"],
])

deploydocs(
    repo="github.com/ImperialCollegeLondon/ReCoDE-NavierStokesPropagator.git",
)
