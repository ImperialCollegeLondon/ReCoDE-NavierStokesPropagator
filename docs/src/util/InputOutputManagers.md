# InputOutputManagers

The InputOutputManagers module is a wrapper around the HDF5.jl I/O library, used for reading and writing results from the simulation.

## Struct
```@docs
InputOutputManager
```

## Public API
```@docs
write_grid(io_m::InputOutputManager, domain::DomainDescriptor{T}) where {T<:AbstractFloat}
write_flowfield(io_m::InputOutputManager, state::State, domain::DomainDescriptor{T}) where {T<:AbstractFloat}
write_flow_statistics(io_m::InputOutputManager, state::State, domain::DomainDescriptor{T}) where {T<:AbstractFloat}
write_attribute(io_m::InputOutputManager, state::State{T}, dt::T, nu::T) where {T<:AbstractFloat}
read_flowfield!(io_m::InputOutputManager, state::State)
```