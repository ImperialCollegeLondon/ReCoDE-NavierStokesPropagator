# SimulationConditions

The SimulationConditions module contain all required information that describes the metadata of the simulation that is unrelated to the domain (DomainDescriptors) or state (State) of the simulation.

## Struct
```@docs
SimulationCondition
```

## Public API
The method parameter are dependent on

    init_cond::Char:
        The initial condition [(p) Poiseuille flow, (c) Couette flow, (r) restart flow from file].

    init_kick::T:
        A floating-point value for the initial kick to the system.

    force_type::Char
        Prescribes the forcing function [(m) mass flow rate, (p) pressure gradient].

    force_magnitude::T
        Prescribes the forcing magnitude in relation to the forcing function.


```@docs
planeCouetteFlow(init_cond::Char, init_kick::T) where {T<:AbstractFloat}
planePoiseuilleFlow(init_cond::Char, init_kick::T, force_type::Char, force_magnitude::T) where {T<:AbstractFloat}
```