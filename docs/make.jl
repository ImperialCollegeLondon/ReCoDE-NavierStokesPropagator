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

makedocs(
    sitename = "Navier-Stokes Propagator",
    pages = [
        "Home" => "index.md",
        "Simulation Constructs" => Any["sim/States.md",],
        "Simulation Utilities" => Any["util/LinearSolvers.md", "util/Transformations.md"],
    ],
)
