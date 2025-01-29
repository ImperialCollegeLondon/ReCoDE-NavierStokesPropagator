# Transformations

The Transformations module includes the definition for a Transformer, and the methods definition for the required discrete Fourier transform, and the fractional and physical grid transformation. 

## Struct
```@docs
Transformer
```

## Public API
The method parameter are dependent on

    tf::Transformer{T}:
        the Transformer struct encapsulating the required FFTW abstractions for the discrete Fourier transform.

    state::State{T}:
        the State{T} struct encapsulating the states of the simulations.

    arr::Array{Complex{T},3}:
        the array to be transformed

    frac_mode::Bool:
        true if the FFTW call is applied on a fractional array

    dealias_mode::Bool:
        true if the 2/3-dealiasing scheme should be applied

```@docs
physical2fourier!(tf::Transformer{T}, state::State{T}) where {T<:AbstractFloat}
physical2fourier!(tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool, dealias_mode::Bool) where {T<:AbstractFloat}
fourier2physical!(tf::Transformer{T}, state::State{T}) where {T<:AbstractFloat}
fourier2physical!(tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool) where {T<:AbstractFloat}
y2y_F!(state::State{T}) where {T<:AbstractFloat}
y_F2y!(state::State)
```