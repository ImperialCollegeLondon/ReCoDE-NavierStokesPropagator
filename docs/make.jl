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

makedocs(sitename = "Navier-Stokes Propagator")
