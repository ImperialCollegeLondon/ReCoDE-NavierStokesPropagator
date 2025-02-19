# Documenter make.jl

#=
Package
=#

include("../src/NavierStokesPropagators.jl")
using .NavierStokesPropagators

using Documenter

#=
Make Documenter
=#

makedocs(;
    sitename="Navier-Stokes Propagator",
    pages=[
        "Home" => Any["index.md", "getting_started.md"],
        "Simulation Constructs" =>
            Any["sim/DomainDescriptors.md", "sim/States.md", "sim/SimulationConditions.md"],
        "Simulation Utilities" => Any[
            "util/Transformations.md",
            "util/LinearSolvers.md",
            "util/InputOutputManagers.md",
        ],
        "Developer Guide" => Any["best_practices.md"],
    ],
)

deploydocs(; repo="github.com/ImperialCollegeLondon/ReCoDE-NavierStokesPropagator.git")
