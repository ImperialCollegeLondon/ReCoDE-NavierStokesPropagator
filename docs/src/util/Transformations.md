# Transformation Documentation

```@meta
CurrentModule = Transformations
```

```@docs
physical2fourier!(tf::Transformer{T}, state::State{T})
```

```@docs
physical2fourier!(tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool, dealias_mode::Bool)
```