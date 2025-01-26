# Transformations

Transformations module defines a Transformer-struct that handles the switching between the physical-fourier domain dimension and the switching of the physical-fractional grids that is defined by the simulation algorithm. It is an utility wrapper around the FFTW package.

```@docs
Transformer
physical2fourier!(tf::Transformer{T}, state::State{T}) where {T<:AbstractFloat}
physical2fourier!(tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool, dealias_mode::Bool) where {T<:AbstractFloat}
fourier2physical!(tf::Transformer{T}, state::State{T}) where {T<:AbstractFloat}
fourier2physical!(tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool) where {T<:AbstractFloat}
y2y_F!(state::State{T}) where {T<:AbstractFloat}
y_F2y!(state::State)
```