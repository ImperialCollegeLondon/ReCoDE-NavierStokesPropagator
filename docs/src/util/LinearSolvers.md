# LinearSolvers

The LinearSolvers module is broken into two components. The first containing the core methods to solve for a TriDiagonal linear system, whilst the second component contain methods to setup the problem specific TriDiagonal system to the simulation to solve for the velocity-pressure divergence updates.

## Struct
```@docs
TriDiagonalMatrix{T<:Number}
LinearSolver{T<:AbstractFloat}
```