# Getting Started

This section details steps to get started with this project, including installation, configuration, and basic usage.

## System

- Install [Julia](https://julialang.org/downloads/) on your system.
- Install [VS Code](https://code.visualstudio.com/Download) on your system.
- Install [Julia language extension](https://code.visualstudio.com/docs/languages/julia) on VS Code.

## Getting Started

1. Clone the repository using `git clone https://github.com/ImperialCollegeLondon/ReCoDE-NavierStokesPropagator.git` and enter the repository directory via `cd ReCoDE-NavierStokesPropagator`.
2. Open the system terminal via VSCode and launch the Julia REPL via `julia`.
3. Launch Julia `pkg` mode via pressing `]` in the Julia REPL.
4. In the `pkg` mode, activate the package environment via `activate .`.
5. Testing can be done by entering `test` in the REPL's pkg mode.

## Quick Start

The first step to running the simulation is to specify the type of simulation configuration to run. Either Plane Couette flow or Plane Poiseulle flow is support in by the package.

### Simulation Configuration

For Plane Poiseulle flow, we have the following configuration:

```julia
# (p) for Plane Poiseulle flow or (r) for restart from file
init_cond::Char = 'p'

# Initial noise perturbation magnitude
init_kick::Float64 = 0.001

# Forcing Constraint - (m) mass flow rate, (p) pressure gradient
force_constraint::Char = 'p'

# Forcing Magnitude
force_magnitude::Float64 = -1.0

sim_cond::SimulationCondition = planePoiseuilleFlow(init_cond, init_kick, force_type, force_magnitude) 
```

For Plane Couette flow, we have the following configuration:

```julia
# (p) for Plane Poiseulle flow or (r) for restart from file
init_cond::Char = 'c'

# Initial noise perturbation magnitude
init_kick::Float64 = 0.001

sim_cond::SimulationCondition = planeCouetteFlow(init_cond, init_kick)
```

### Simulation Domain

For the simulation domain, the configuration is set by:

```julia
# Domain Size
Lx::Float64, Ly::Float64, Lz::Float64 = (1.75pi, 2.0, 1.2pi)

# Mesh Grid Size
(nx::Int64, ny::Int64, nz::Int64) = (33, 65, 33)

domain::DomainDescriptor = DomainDescriptor{Float64}(nx, ny, nz, Lx, Ly, Lz)
```

### Temporal Parameter

For the temporal parameter, we created a *time-keeper* struct that keeps timesteps in sync.

```julia
# Total Simulation Time
T_sim::Float64 = 5000.0

# Time Step
dt::Float64 = 0.05

# Time Step Bound
dt_max::Float64 = 0.09
dt_min::Float64 = 0.0005

# Maximum No. of Time Step
nt_max::Int64 = floor(Int64, T_sim / dt_min) + 1

# CFL
cfl::Float64 = 0.4

# Adaptive Timestep
adaptive_dt::Bool = false

# Time-Keeper
t_stepper::TimeStepper = TimeStepper{Float64}(dt, adaptive_dt, dt_max, dt_min, nt_max, cfl)

# State
state::State = State{FloatType}(T_sim, nx, ny, nz)
```

### Utility Component
The 2 utility components, `InptOutputManager` and `Transformer`, are initailised via:

```julia
# State I/O (n time-step)
freq_state::Int64 = 1000

# Statistics I/O (n time-step)
freq_stats::Int64 = 200
min_step_stats::Int64 = 1000

# I/O Manager
io_m::InputOutputManager = InputOutputManager("input", "output", freq_state, freq_stats, min_step_stats)

# Transformer
tf::Transformer = Transformer{FloatType}("input/fftw_wisdom", state, domain; wisdom_flag="exhaustive")
```

### Run Simulation
Lastly, we can run the simulation via:

```julia
ns::NavierStokesPropagator = NavierStokesPropagator(t_stepper, domain, state, sim_cond, io_m, tf, FloatType)

run_simulation!(ns)
```

## FFTW Abstraction

In the core simulation algorithm, the velocity-pressure states are required to be transformed between physical and fourier domain. This is where the FFTW (Fast Fourier Transform in the West) comes into play. Since the transformation is done in-place for each array, the main way the FFTW library interfaces with these arrays are with a FFTW plan.

```julia
fft_plan * array
```

However, a distinction is required to be made for the discrete Fourier coefficients with the Fourier series coefficient which has a scaling of `nx * nz` between them. So additionally to a transformation operation as in the code snippet above, a normalisation of `nx * nz` is required as well.

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

The combination of all these optional and package specific features leads to the public APIs written for the `Transformations` module, allowing the develop to repeatably use the following transformation function as ease.

```julia
physical2fourier!(
    tf::Transformer{T}, arr::Array{Complex{T},3}, frac_mode::Bool, dealias_mode::Bool
)
```

## HDF5 Abstraction

For input/output of velocity-pressure state arrays, the HDF5 package is used with its defined groups and dataset formalism. The public APIs written for the `InputOutputManagers` module allows the group definition and datasets to be written out in code, having a clear separation of responsibilities, segregated by modules. Hence, for a user reading code in `run_simulation!(ns::NavierStokesPropagator)`, `write_flowfield(io_m::InputOutputManager)` would clearly been seen as an abstraction of output. The alternative way of not writing the `InputOutputManagers` wrapper module would required simulation code and HDF5 API calls to be interleaved in `run_simulation!(ns::NavierStokesPropagator)`, leading to convoluted and harder-to-read code.

A key point to make here is the switching of the `x` and `y` dimensions in the simulation code for 3D arrays. In the simulation, since we predominantly use array operators in `y` and loop through `x` and `z` dimensions, it was sensible for arrays to be arranged in `(y,x,z)`, allowing for contiguous memory layout to be achieved. However, normal convention would expect for a `(x,y,z)` layout, hence, this extra layer of complication can be readily address by the `InputOutputManagers` module which is responsible for the HDF5 abstraction. As seen in `read_flowfield!(io_m::InputOutputManger, state::State)`, we have calls to `PermutedDimsArray` handling the permuation of dimensions for the user.

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