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
# Kinematic Viscosity
nu::Float64 = 1.0/2200.0

# (p) for Plane Poiseulle flow or (r) for restart from file
init_cond::Char = 'p'

# Initial noise perturbation magnitude
init_kick::Float64 = 0.001

# Forcing Constraint - (m) mass flow rate, (p) pressure gradient
force_constraint::Char = 'p'

# Forcing Magnitude
force_magnitude::Float64 = -1.0

sim_cond::SimulationCondition = planePoiseuilleFlow(init_cond, init_kick, force_constraint, force_magnitude, nu) 
```

For Plane Couette flow, we have the following configuration:

```julia
# Kinematic Viscosity
nu::Float64 = 1.0/400.0

# (c) for Plane Couette flow or (r) for restart from file
init_cond::Char = 'c'

# Initial noise perturbation magnitude
init_kick::Float64 = 0.001

sim_cond::SimulationCondition = planeCouetteFlow(init_cond, init_kick, nu)
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
state::State = State{Float64}(T_sim, nx, ny, nz)
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
tf::Transformer = Transformer{Float64}("input/fftw_wisdom", state, domain; wisdom_flag="exhaustive")
```

### Run Simulation
Lastly, we can run the simulation via:

```julia
ns::NavierStokesPropagator = NavierStokesPropagator(t_stepper, domain, state, sim_cond, io_m, tf, Float64)

run_simulation!(ns)
```