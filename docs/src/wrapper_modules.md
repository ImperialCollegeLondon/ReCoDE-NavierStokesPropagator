# Wrapper Modules

## FFTW Abstraction

In the core simulation algorithm, the velocity-pressure states are required to be transformed between physical and fourier domain. This is where the Fast Fourier Transform in the West (FFTW) comes into play. Since the transformation is done in-place for each array, the main way the FFTW library interfaces with these arrays are with a FFTW plan.

```julia
fft_plan * array
```

However, a distinction is required to be made for the discrete Fourier coefficients with the Fourier series coefficient which has a scaling of `nx * nz` between them. So in addition to a transformation operation as shown in the code snippet above, a normalisation of `nx * nz` is required as well.

```julia
fft_plan * array

array[:] ./= nx * nz
```
For the treatment of non-linear terms, an additional dealiasing-mode is required, leading to the zero-ing of frequencies larger than a threshold.

```julia
# De-aliasing
if (dealias_mode)
    for ix in 1:nx
        if (abs(tf.freq_kx[ix]) > tf.freq_kx_max)
            arr[:, ix, :] .= 0.0 + 0.0im
        end
    end

    for iz in 1:nz
        if (abs(tf.freq_kz[iz]) > tf.freq_kz_max)
            arr[:, :, iz] .= 0.0 + 0.0im
        end
    end
end
```

The combination of all these optional and package specific features leads to the public APIs written for the `Transformations` module, allowing the developer to repeatably use the following transformation function as ease.

```julia
physical2fourier!(
    tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool, dealias_mode::Bool
)
```

## HDF5 Abstraction

For input/output of velocity-pressure state arrays, the HDF5 package is used with its defined groups and dataset formalism. The public APIs written for the `InputOutputManagers` module allows the group definition and datasets to be written out in code, having a clear separation of responsibilities, segregated by modules. Hence, for a user reading code in `run_simulation!(ns::NavierStokesPropagator)`, `write_flowfield(io_m::InputOutputManager)` would clearly been seen as an abstraction of output. The alternative way of not writing the `InputOutputManagers` wrapper module would require simulation code and HDF5 API calls to be interleaved in `run_simulation!(ns::NavierStokesPropagator)`, leading to convoluted and harder-to-read code.

A key point to make here is the switching of the `x` and `y` dimensions in the simulation code for 3D arrays. In the simulation, since we predominantly use array operators in `y` and loop through `x` and `z` dimensions, it makes sense for the arrays to be arranged in `(y,x,z)`, allowing for contiguous memory layout to be achieved. However, normal convention would expect for a `(x,y,z)` layout, hence, this extra layer of complication can be readily addressed by the `InputOutputManagers` module which is responsible for the HDF5 abstraction. As seen in `read_flowfield!(io_m::InputOutputManger, state::State)`, we have calls to `PermutedDimsArray` handling the permutation of dimensions for the user.

```julia
# Read State and Attribute
h5open(input_path, "r") do fid
    g::HDF5.Group = open_group(fid, "/")

    # U Velocity
    state.u[:, :, :] = PermutedDimsArray(read(g, "u"), (2, 1, 3))

    # V Velocity
    state.v[:, :, :] = PermutedDimsArray(read(g, "v"), (2, 1, 3))

    # W Velocity
    state.w[:, :, :] = PermutedDimsArray(read(g, "w"), (2, 1, 3))
end
```