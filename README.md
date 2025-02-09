# Navier-Stokes Propagator

## Description

The repository contains source files for a Julia package written to solve the chaotic Navier-Stokes equations. The project aims to exemplify the improvement of simulation code quality using good programming abstractions (Julia Structs) and defining sensible simulation constructs that reflect each core component of the simulation. Additionally, the repository adopts good programming practices by implementing detailed documentation with tests written for each simulation construct.

The algorithm used here is inspired from [Diablo](https://github.com/johnryantaylor/DIABLO), written by John R. Taylor, with minor modifications to the mesh spacing.

## Learning Outcomes

- Ability to distinguish and compartmentalise a multi-facet simulation algorithm into sensible simulation constructs via the implementation with Julia Structs 
- Adopting good testing methodologies, through implementation of unit and integration tests
- Understanding the FFTW and HDF5 libraries and their abstraction through the usage of utility modules in their respective source files

| Task       | Time    |
| ---------- | ------- |
| Reading    | 3 hours |

## Requirements

### Theoretical
- Exploring Fourier series and its relation to the discrete Fourier transform
- Understanding of the finite difference method for approximating spatial derivatives and the use of implicit/explicit time-stepping schemes
- Familiarity with Julia. The Early Career Researcher Institute also provides an [Introduction to Julia](https://www.imperial.ac.uk/students/academic-support/graduate-school/professional-development/doctoral-students/research-computing-data-science/courses/introduction-to-julia/) course.

Note : The theoretical understanding of the finite difference method, Fourier series and the time-stepping scheme are not a major necessity. An understanding of the mathematical structure for algorithimic implementation is sufficient to understand the ideas discussed in this tutorial.

### System

- Install [Julia](https://julialang.org/downloads/) on your system.
- Install [VS Code](https://code.visualstudio.com/Download) on your system.
- Install [Julia language extension](https://code.visualstudio.com/docs/languages/julia) on VS Code.

## Getting Started

1. Clone the repository using `git clone https://github.com/ImperialCollegeLondon/ReCoDE-NavierStokesPropagator.git` and enter the repository directory via `cd ReCoDE-NavierStokesPropagator`.
2. Open the system terminal via VSCode and launch the Julia REPL via `julia`.
3. Launch Julia `pkg` mode via pressing `]` in the Julia REPL.
4. In the `pkg` mode, activate the package environment via `activate .`.
5. Testing can be done by entering `test` in the REPL's pkg mode.

## Project Organisation

Solving the Navier-Stokes equations with periodic boundary conditions in the streamwise and spanwise directions, with walls bounding the wall-normal direction, involves the use of the Fourier-Galerkin technique. That is to say that the numerics involved in the time-stepping algorithm will transform between the physical and frequency domain during the numerical solve for the turbulent boundary layer. 

From the perspective of software, we can break the neatly into two components. The first are termed `Simulation Constructs`, essentially the constructs of the program that store the direct physical data / variables involved during the simulation. On the other hand, `Simulation Util` are defined to be constructions of the program that aids with the simulation that are not directly related to the physical problem at hand. Normally, these source files consist of wrapper function that wraps around open-source utility libraries to further simplify their utility, bringing these libraries closer towards a more user friendly experience Being guided by the principle of minimising the coupling between these source files, `Simulation Constructs` and `Simulation Util` are defined in the following project structure (see below). In other words, the followng set of division of source code reduces the cross-dependencies that are required during the design of the codebase.

`Simulation Constructs` acts as the foundation of the codebase, with the two primary Julia-structs, `DomainDescriptors.jl` and `States.jl`, being stand-alone concepts. DomainDescriptors contain all required mesh information and the sizing of the computational box, while States itself encapsulate the velocity fields and the state of simulation that are mutable, and where the progression of the time-stepping scheme is dependent on. In a sense, these constructs are the conceptually the closest to the *physical* problem at hand. Building on top of these two constructs, `SimulationConditions.jl` encapsulates the initialisation of the simulation and the boundary conditions the simulation is required to adhere to. Although `SimulationConditions.jl` depend on `DomainDescriptors.jl` and `States.jl`, the following conceptual division leads to a minimisation of the dependencies required. Looking at the code base, the initialisation of each Julia-struct can be done by referring to other variables/functions within its initialiser.

### Project Structure
```log
.
├── src
│   ├── NavierStokesPropagators.jl
│   ├── NavierStokesPropagatorsCore.jl
│   ├── sim
│   │   ├── DomainDescriptors.jl
│   │   ├── SimulationConditions.jl
│   │   └── States.jl
│   └── util
│       ├── InputOutputManagers.jl
│       ├── LinearSolvers.jl
│       └── Transformations.jl
├── test
│   ├── sim (mirrors src/sim)
│   ├── util (mirrors src/util)
│   ├── StateTransformations.jl
│   └── runtest.jl
├── docs
└── .github
    └── (.yml files for CI/CD)
```

## License

This project is licensed under the BSD-3-Clause license.
