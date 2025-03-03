# Best Practices

This section details the software engineering best practices adopted for this project, including automated linting and formatting, automated unit testing, and automated publishing of dynamic/static documentation - all using native Julia libraries!

## Project Organisation

Solving the Navier-Stokes equations with periodic boundary conditions in the streamwise and spanwise directions, with walls bounding the wall-normal direction, involves the use of the Fourier-Galerkin technique. That is to say that the numerics involved in the time-stepping algorithm will transform between the physical and frequency domain during the numerical solve for the turbulent boundary layer. 

From the perspective of software, we can break the neatly into two components. The first are termed `Simulation Constructs`, essentially the constructs of the program that store the direct physical data / variables involved during the simulation. On the other hand, `Simulation Util` are defined to be constructions of the program that aids with the simulation that are not directly related to the physical problem at hand. Normally, these source files consist of wrapper function that wraps around open-source utility libraries to further simplify their utility, bringing these libraries closer towards a more user friendly experience Being guided by the principle of minimising the coupling between these source files, `Simulation Constructs` and `Simulation Util` are defined in the following project structure (see below). In other words, the followng set of division of source code reduces the cross-dependencies that are required during the design of the codebase.

![Dependency Level](../assets/dep.png "Dependency Hierachy")

`Simulation Constructs` acts as the foundation of the codebase, with the two primary Julia-structs, `DomainDescriptors.jl` and `States.jl`, being stand-alone concepts. DomainDescriptors contain all required mesh information and the sizing of the computational box, while States itself encapsulate the velocity fields and the state of simulation that are mutable, and where the progression of the time-stepping scheme is dependent on. In a sense, these constructs are the conceptually the closest to the *physical* problem at hand. Building on top of these two constructs, `SimulationConditions.jl` encapsulates the initialisation of the simulation and the boundary conditions the simulation is required to adhere to. Looking at the code base, the initialisation of each Julia-struct can be done by referring to other variables/functions within its initialiser.

`Simulation Util` instead, as discussed above, wraps around open-source libraries and algorithmic routines that aids in the simulation of the Navier-Stokes equations. `Simulation Constructs` and `Simulation Util` are divided in a way that source files in `Simulation Util` would only be dependent on `Simulation Constructs` but not vis versa.

Dependency of InputOutputManager.jl and Transformations.jl
```julia
using ..States: State
using ..DomainDescriptors: DomainDescriptor
```

## Testing

Here, we introduce the testing methodology for software, based on the `Test.jl` package in Julia's standard library. As software gets developed, modifications to the codebase would inevitably change the behaviour of outputs. By placing in safeguards in terms of automated testing (since it is software that can be run automatically by a microprocessor), bugs introduced during the development can be avoided and flagged to the developer in a manner appropriate.

In general, testing can be broken into levels of granularity. At the lowest level of granularity, we have unit tests, suitably written for testing the small chunks of code, generally functions and basic properties of structures. Climbing one level up, we have integration testing. Here, the focus is more on validation of correct software behaviour between inputs and outputs of functions. As the software written gets more complicated, systems level testing is introduced as the next level of granularity, which now is considered a black-box technique, usually tasked to be executed by software quality engineers instead of developers and hence we omit ourselves from the systems testing and higher level granularity at this point.

In terms of the unit tests written, coverage has been written for `DomainDesriptors.jl`, `States.jl`, `LinearSolvers.jl` and `Transformations.jl`. Due to its nature of being a unit test, this package adopts the convention that the test files are named similarly to their respective source files, and structure in a mirrored way to `src`. In Julia's `Test.jl` package, developers are able to assign test cases and apply assertions via the `@test` macro.
```julia
@test <write statement here>
```
This is the core feature of testing, allowing for assertion to be made to ensure the correctness of the statement in the <> braces. Doing so in `test/sim/DomainDescriptors.jl`, we are able to test the constructor of `DomainDescriptor`, ensuring the features expected from the `DomainDescriptor` struct are adhered to.

Furthermore, due to the multiple assertions being made in a test, the `@testset` macro allows the developer to sort sets together, organising the test summary into digestable chunks. When organised wisely where each testset are sorted by common features, each testset can elucidate insight for the starting point of debugging the codebase when bugs are detected from testing.

Going one level up, the testing done by `StateTransformations.jl` starts integrating the properties of `State` which are to be transformed by properties in a `Transformer`. Since the interaction here are due to functions involving a `util` and a `sim` dependency, it is considered as an integration test, and have resides in the `StateTransformations.jl` test file in `test/`.

### Unit Testing Code Coverage

Unit tests should ideally be written for all cases, having coverage for all lines of the source code. However, it is impractical as exhibit in this instance. The reader may notice that although the unit tests have been written for most of the `Simulation Construct (\sim)` and `Simulation Utility (\util)` source files, none have been written for the `NavierStokesPropagators*.jl` source files. The reason is the unavailability of good quality tests to be written as the simulation code itself are complicated and involved. We are dealing with the numerics of chaotic turbulence here! Hence, as a developer, it is important to realise the issue apriori and include safeguards into the codebase to guard from introducing bugs.

For this simulation, the state arrays are be in the `frac_mode` and `freq_mode` (essentially in the frequency domain). However, the hundreds of lines of code obfuscate the currently state (`frac_mode` and `freq_mode`) when a function in `NavierStokesPropagators*.jl` runs. As a result, checks are placed as safeguards to the entry of a function, for example,
```julia
# State Mode Check
if (state.fourier_mode == false)
    @error("Computing flow attribute when fourier_mode = $(state.fourier_mode).")

    return nothing
end
```
to avoid developers from entering the function in the wrong states.
