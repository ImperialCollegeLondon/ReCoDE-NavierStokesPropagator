<!--
This includes your top-level README as you index page i.e. homepage.

This will not be the best approach for all exemplars, so feel free to customise
your index page as you see fit.
-->

{%
include-markdown "../README.md"

%}

```@docs
physical2fourier!(tf::Transformer{T}, state::State{T}) where {T<:AbstractFloat}
```

```@docs
physical2fourier!(tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool, dealias_mode::Bool) where {T<:AbstractFloat}
```

<!-- Add more files in the `docs/` directory for them to be automatically
included in the Mkdocs documentation -->
