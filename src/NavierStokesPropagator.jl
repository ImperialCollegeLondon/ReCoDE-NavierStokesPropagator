module NavierStokesPropagator

#=
Package
=#

include("sim/DomainDescriptors.jl")
using .DomainDescriptors

include("sim/States.jl")
using .States

include("sim/SimulationConditions.jl")
using .SimulationConditions

include("util/InputOutputManagers.jl")
using .InputOutputManagers

include("util/Transformations.jl")
using .Transformations

include("util/LinearSolvers.jl")
using .LinearSolvers

#=
Public API
=#

export BoundTuple
export SimulationCondition, planeCouetteFlow, planePoiseuilleFlow

export Transformer
export y2y_F!, y_F2y!
export physical2fourier!, fourier2physical!

export TriDiagonalMatrix, LinearSolver
export solve_system_R!, solve_system_C!

export InputOutputManager
export write_grid
export write_flowfield, read_flowfield!
export write_flow_statistics, write_attribute

export NavierStokesPropagator, TimeStepper
export DomainDescriptor, SimulationCondition, State
export InputOutputManager, Transformer

export run_simulation!

#=
Navier Stokes Propagator
=#

mutable struct TimeStepper{T<:AbstractFloat}

    # Stepsize
    dt::T

    # Adaptive Stepsize
    bool_adaptive::Bool
    dt_max::T # Max Stepsize
    dt_min::T # Min Stepsize
    nt_max::Int64 # Max No. of Timestep
    cfl::T # CFL Condition

    # Runge-Kutta-Wray Scheme
    h_bar::Vector{T}
    beta::Vector{T}
    zeta::Vector{T}

    function TimeStepper{T}(dt::T, bool_adaptive::Bool, dt_max::T, dt_min::T, nt_max::Int64, cfl::T) where {T<:AbstractFloat}

        h_bar::Vector{T} = [8.0 / 15.0; 2.0 / 15.0; 5.0 / 15.0]
        beta::Vector{T} = [1.0; 25.0 / 8.0; 9.0 / 4.0]
        zeta::Vector{T} = [0.0; -17.0 / 8.0; -5.0 / 4.0]

        return new(dt, bool_adaptive, dt_max, dt_min, nt_max, cfl, h_bar, beta, zeta)

    end

end

struct NavierStokesPropagator

    # Simulation Construct
    t_stepper::TimeStepper
    domain::DomainDescriptor
    state::State
    sim_cond::SimulationCondition

    # Utility Construct
    solver::LinearSolver
    io_m::InputOutputManager
    tf::Transformer

    # FloatType
    FloatType::Type

    function NavierStokesPropagator(t_stepper::TimeStepper, domain::DomainDescriptor, state::State, sim_cond::SimulationCondition,
        io_m::InputOutputManager, tf::Transformer, FloatType::Type)

        solver::LinearSolver{FloatType} = LinearSolver{FloatType}(domain.ny, domain.ny_F)

        return new(t_stepper, domain, state, sim_cond, solver, io_m, tf, FloatType)

    end

end

#=
Utility Algorithm
=#

# Flow Attribute Computation
function update_flow_attribute!(ns::NavierStokesPropagator)

    # Simulation Construct
    domain::DomainDescriptor = ns.domain
    state::State = ns.state
    sim_cond::SimulationCondition = ns.sim_cond

    # State Mode Check
    if (state.fourier_mode == false)

        @error("Computing flow attribute when fourier_mode = $(state.fourier_mode).")

        return

    end

    # Bulk Velocity and Shear Velocity
    state.u_bulk = integral_y(state.u[:, 1, 1], domain.dy) / (domain.y[end] - domain.y[1])

    tau_lwr::Float64 = sim_cond.nu * (real(state.u[2, 1, 1] - state.u[1, 1, 1])) / domain.dy[1]
    tau_upr::Float64 = -sim_cond.nu * (real(state.u[end, 1, 1] - state.u[end-1, 1, 1])) / domain.dy[end]

    # Friction Velocity
    state.u_tau_lwr = sqrt(abs(tau_lwr))
    state.u_tau_upr = sqrt(abs(tau_upr))

    # Pressure Gradient Update
    if (sim_cond.force_constraint == 'm')

        state.dp0_dx = -(tau_lwr + tau_upr) / domain.Ly
        state.dp0_dx += 1.5 * (state.u_bulk - sim_cond.force_magnitude)

    end

    return nothing

end

# Flow Statistics Computation
function update_flow_statistics!(ns::NavierStokesPropagator)

    # Simulation Construct
    domain::DomainDescriptor = ns.domain
    state::State = ns.state

    # State Mode Check
    if (state.fourier_mode == false)

        @error("Updating flow statistics when fourier_mode = $(state.fourier_mode).")

        return

    end

    # Update Sample Size
    state.n_sample += 1

    # Temporal Mean Velocity
    state.u_bar[:] += real(state.u[:, 1, 1] - state.u_bar[:]) / state.n_sample
    state.v_bar[:] += real(state.v[:, 1, 1] - state.v_bar[:]) / state.n_sample
    state.w_bar[:] += real(state.w[:, 1, 1] - state.w_bar[:]) / state.n_sample

    # Transformation Operation
    fourier2physical!(ns.tf, ns.state)

    # Temporal Mean Reynolds Stress
    uu_bar::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)
    vv_bar::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)
    ww_bar::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)

    uv_bar::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)
    uw_bar::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)
    vw_bar::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)

    for iz in eachindex(domain.z)
        for ix in eachindex(domain.x)

            uu_bar[:] += (state.u[:, ix, iz] .- state.u_bar[:]) .* (state.u[:, ix, iz] .- state.u_bar[:])
            vv_bar[:] += (state.v[:, ix, iz] .- state.v_bar[:]) .* (state.v[:, ix, iz] .- state.v_bar[:])
            ww_bar[:] += (state.w[:, ix, iz] .- state.w_bar[:]) .* (state.w[:, ix, iz] .- state.w_bar[:])

            uv_bar[:] += (state.u[:, ix, iz] .- state.u_bar[:]) .* (state.v[:, ix, iz] .- state.v_bar[:])
            uw_bar[:] += (state.u[:, ix, iz] .- state.u_bar[:]) .* (state.w[:, ix, iz] .- state.w_bar[:])
            vw_bar[:] += (state.v[:, ix, iz] .- state.v_bar[:]) .* (state.w[:, ix, iz] .- state.w_bar[:])

        end
    end

    state.uu_bar[:] += (uu_bar[:] / (domain.nx * domain.nz) - state.uu_bar[:]) / state.n_sample
    state.vv_bar[:] += (vv_bar[:] / (domain.nx * domain.nz) - state.vv_bar[:]) / state.n_sample
    state.ww_bar[:] += (ww_bar[:] / (domain.nx * domain.nz) - state.ww_bar[:]) / state.n_sample

    state.uv_bar[:] += (uv_bar[:] / (domain.nx * domain.nz) - state.uv_bar[:]) / state.n_sample
    state.uw_bar[:] += (uw_bar[:] / (domain.nx * domain.nz) - state.uw_bar[:]) / state.n_sample
    state.vw_bar[:] += (vw_bar[:] / (domain.nx * domain.nz) - state.vw_bar[:]) / state.n_sample


    # Transformation Operation
    physical2fourier!(ns.tf, ns.state)

    return nothing

end

# CFL Time-step Update
function timestep_refinement!(ns::NavierStokesPropagator)

    @info("Computing timestep refinement.")

    # Simulation Construct
    t_stepper::TimeStepper = ns.t_stepper
    domain::DomainDescriptor = ns.domain
    state::State = ns.state

    # Mesh Dependent Step Size
    # ds::ns.FloatType = min(domain.dx, domain.dy[1])
    # dt_max::ns.FloatType = 0.5 * ds^2 / ns.sim_cond.nu

    # Velocity Dependent Step Size
    dt::ns.FloatType = 999.0

    dt_x::ns.FloatType = 999.0
    dt_y::ns.FloatType = 999.0
    dt_z::ns.FloatType = 999.0

    for iz in eachindex(domain.z)
        for ix in eachindex(domain.x)
            for iy = 3:domain.ny-2

                # Time Step
                if (abs(state.u[iy, ix, iz]) > 1e-3)
                    dt_x = t_stepper.cfl * domain.dx / abs(state.u[iy, ix, iz])
                end

                if (abs(state.v[iy, ix, iz]) > 1e-3)
                    dt_y = t_stepper.cfl * domain.dy_F[iy] / abs(state.v[iy, ix, iz])
                end

                if (abs(state.w[iy, ix, iz]) > 1e-3)
                    dt_z = t_stepper.cfl * domain.dz / abs(state.w[iy, ix, iz])
                end

                if (min(dt_x, dt_y, dt_z) < dt)
                    dt = min(dt_x, dt_y, dt_z)
                end

            end
        end
    end

    # Time Stepper
    dt_prev::ns.FloatType = t_stepper.dt

    t_stepper.dt = dt

    # if (dt_max > dt) && (ns.sim_cond.init_cond == 'p')
    #     t_stepper.dt = dt_max
    # else
    #     t_stepper.dt = dt
    # end

    # Bounding
    if (t_stepper.dt > t_stepper.dt_max)
        t_stepper.dt = t_stepper.dt_max
    end

    if (t_stepper.dt < t_stepper.dt_min)
        t_stepper.dt = t_stepper.dt_min
    end

    @info("Timestep refinement applied: dt = $(dt_prev) -> $(t_stepper.dt)")

    return nothing

end

#=
Initial Conditions
=#

function init_flowfield!(sim_cond::SimulationCondition{T}, domain::DomainDescriptor{T}, state::State{T}, tf::Transformer{T}) where {T<:AbstractFloat}

    # Apply Perturbation
    apply_perturbation!(sim_cond, domain, state, tf)

    # Plane Poiseuille Flow
    if (sim_cond.init_cond == 'p')

        if (sim_cond.force_constraint == 'm')

            for iy = 1:domain.ny
                state.u[iy, 1, 1] += (1.0 - domain.y[iy]^2) * sim_cond.force_magnitude / (2.0 / 3.0)
            end

        elseif (sim_cond.force_constraint == 'p')

            for iy = 1:domain.ny
                state.u[iy, 1, 1] += (1.0 / sim_cond.nu) * (1.0 - domain.y[iy]^2) * (-0.5 * sim_cond.force_magnitude)
            end

        end

    end

    # Plane Couette Flow
    if (sim_cond.init_cond == 'c')

        for iy = 1:domain.ny
            state.u[iy, 1, 1] += domain.y[iy]
        end

    end

    # Enforce Boundary Condition
    state.u[1, :, :] .= 0.0
    state.u[1, 1, 1] = sim_cond.bound_cond[(vel='u', wall='l')]

    state.u[end, :, :] .= 0.0
    state.u[end, 1, 1] = sim_cond.bound_cond[(vel='u', wall='u')]

    state.v[1, :, :] .= 0.0
    state.v[1, 1, 1] = sim_cond.bound_cond[(vel='v', wall='l')]

    state.v[end, :, :] .= 0.0
    state.v[end, 1, 1] = sim_cond.bound_cond[(vel='v', wall='u')]

    state.w[1, :, :] .= 0.0
    state.w[1, 1, 1] = sim_cond.bound_cond[(vel='w', wall='l')]

    state.w[end, :, :] .= 0.0
    state.w[end, 1, 1] = sim_cond.bound_cond[(vel='w', wall='u')]

    return nothing

end

function apply_perturbation!(sim_cond::SimulationCondition{T}, domain::DomainDescriptor{T}, state::State{T}, tf::Transformer{T}) where {T<:AbstractFloat}

    rng = Xoshiro(Int64(floor(time())))

    # Apply Perturbation
    iz_c::Int64 = 0
    ix_c::Int64 = 0

    for iz = 2:domain.nz

        freq_kz::T = tf.freq_kz[iz]

        for fz = 1:domain.nz

            if (tf.freq_kz[fz] == -freq_kz)
                iz_c = fz
            end

        end

        for ix = 1:domain.nx

            freq_kx::T = tf.freq_kx[ix]

            for fx = 1:domain.nx

                if (tf.freq_kx[fx] == -freq_kx)
                    ix_c = fx
                end

            end

            for iy = 2:domain.ny-1

                u::T = (rand(rng, T) - 0.5) * sim_cond.init_kick / 2.0
                state.u[iy, ix, iz] += u
                state.u[iy, ix_c, iz_c] += u

                v::T = (rand(rng, T) - 0.5) * sim_cond.init_kick / 2.0
                state.v[iy, ix, iz] += v
                state.v[iy, ix_c, iz_c] += v

                w::T = (rand(rng, T) - 0.5) * sim_cond.init_kick / 2.0
                state.w[iy, ix, iz] += w
                state.w[iy, ix_c, iz_c] += w

            end

        end

    end

    return nothing

end

function integral_y(arr::Vector{Complex{T}}, dy::Vector{T})::Complex{T} where {T<:AbstractFloat}

    # Integral
    integral::Complex{T} = 0.0

    # Loop
    for i = 1:length(dy)

        integral += dy[i] * (arr[i+1] + arr[i]) / 2.0

    end

    return real(integral)

end

#=
Core Algorithm
=#

include("NavierStokesPropagatorsCore.jl")

end
