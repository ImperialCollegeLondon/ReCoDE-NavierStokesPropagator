# Core Time Stepping Algorithm

#=
Core Algorithm
=#

# Dynamics Vector
struct DynamicsArray{T<:AbstractFloat}

    # Runge-Kutta Terms
    RK_U::Array{Complex{T},3}
    RK_V::Array{Complex{T},3}
    RK_W::Array{Complex{T},3}

    # Crank-Nicolson Terms
    IMEX_U::Array{Complex{T},3}
    IMEX_V::Array{Complex{T},3}
    IMEX_W::Array{Complex{T},3}

    # Non-linear Advection
    S::Array{Complex{T},3}
    S_F::Array{Complex{T},3}

    function DynamicsArray{T}(domain::DomainDescriptor) where {T<:AbstractFloat}

        # Grid Size
        nx::Int64 = domain.nx
        nz::Int64 = domain.nz
        ny::Int64 = domain.ny
        ny_F::Int64 = domain.ny_F

        grid_base::NTuple{3,Int64} = (ny, nx, nz)
        grid_frac::NTuple{3,Int64} = (ny_F, nx, nz)

        # Allocation
        RK_U::Array{Complex{T},3} = zeros(Complex{T}, grid_base)
        RK_V::Array{Complex{T},3} = zeros(Complex{T}, grid_frac)
        RK_W::Array{Complex{T},3} = zeros(Complex{T}, grid_base)

        IMEX_U::Array{Complex{T},3} = zeros(Complex{T}, grid_base)
        IMEX_V::Array{Complex{T},3} = zeros(Complex{T}, grid_frac)
        IMEX_W::Array{Complex{T},3} = zeros(Complex{T}, grid_base)

        S::Array{Complex{T},3} = zeros(Complex{T}, grid_base)
        S_F::Array{Complex{T},3} = zeros(Complex{T}, grid_frac)

        return new(RK_U, RK_V, RK_W, IMEX_U, IMEX_V, IMEX_W, S, S_F)

    end

end

# Run Simulation
function run_simulation!(ns::NavierStokesPropagator)

    # Dynamics Array
    R::DynamicsArray = DynamicsArray{ns.FloatType}(ns.domain)

    # Initialisise Initial Condition
    init_simulation!(ns)

    # Output Grid
    write_grid(ns.io_m, ns.domain)

    # Transformation Operation
    fourier2physical!(ns.tf, ns.state)

    # State Initial Condition I/O
    y_F2y!(ns.state)
    write_flowfield(ns.io_m, ns.state, ns.domain)
    y2y_F!(ns.state)

    # Transformation Operation
    physical2fourier!(ns.tf, ns.state)

    # Flow Attribute I/O
    update_flow_attribute!(ns)
    write_attribute(ns.io_m, ns.state, ns.t_stepper.dt, ns.sim_cond.nu)

    # Time-stepping
    while (ns.state.t_sim < ns.state.T_sim) && (ns.state.nt < ns.t_stepper.nt_max)

        # Timestepping Computation
        step!(ns, R)

        # Update Time Variable
        ns.state.nt += 1
        ns.state.t_sim = ns.t_stepper.dt * ns.state.nt

        #=
        I/O Operation
        =#

        # Flow Attribute I/O
        update_flow_attribute!(ns)

        if (mod(ns.state.nt, 5) == 0)

            write_attribute(ns.io_m, ns.state, ns.t_stepper.dt, ns.sim_cond.nu)
            flush(log_io)

        end

        # State / Statistics I/O
        bool_stats_io::Bool =
            (mod(ns.state.nt, ns.io_m.freq_stats) == 0) &&
            (ns.state.nt > ns.io_m.min_step_stats)
        bool_state_io::Bool = (mod(ns.state.nt, ns.io_m.freq_state) == 0)

        if (bool_state_io) || (bool_stats_io)

            y_F2y!(ns.state)

            if (bool_state_io)

                # Transformation Operation
                fourier2physical!(ns.tf, ns.state)

                # Output Flowfield
                write_flowfield(ns.io_m, ns.state, ns.domain)

                # Transformation Operation
                physical2fourier!(ns.tf, ns.state)

            end

            if (bool_stats_io)

                # Time Refinement
                if (ns.t_stepper.bool_adaptive)

                    timestep_refinement!(ns)

                end

                # Output Statistics
                update_flow_statistics!(ns)
                write_flow_statistics(ns.io_m, ns.state, ns.domain)

            end

            y2y_F!(ns.state)

        end

    end

    return nothing

end

# Simulation Initialisation
function init_simulation!(ns::NavierStokesPropagator)

    # Initalise States
    if (ns.sim_cond.init_cond == 'r')

        read_flowfield!(ns.io_m, ns.state)
        physical2fourier!(ns.tf, ns.state)

    else

        physical2fourier!(ns.tf, ns.state)

        init_flowfield!(ns.sim_cond, ns.domain, ns.state, ns.tf)

        fourier2physical!(ns.tf, ns.state)
        physical2fourier!(ns.tf, ns.state)

    end

    # Transformation Operation
    y2y_F!(ns.state)

    # Pressure Update
    continuity_enforcement!(ns, 0)

    # Pressure Poisson Solve
    pressure_poisson_solve(ns)

    # Set Constant Pressure Gradient
    if (ns.sim_cond.force_constraint == 'p')

        ns.state.dp0_dx = ns.sim_cond.force_magnitude

    end

    return nothing

end

# Timestepping Computation
function step!(ns::NavierStokesPropagator, R::DynamicsArray)

    # Simulation Construct
    t_stepper::TimeStepper = ns.t_stepper
    domain::DomainDescriptor = ns.domain
    sim_cond::SimulationCondition = ns.sim_cond
    state::State = ns.state

    # Multi-Stage Time Stepping
    dt::ns.FloatType = t_stepper.dt

    # Wall-Normal Array Bound
    ny::Int64 = domain.ny
    ny_F::Int64 = domain.ny_F

    for i_rk = 1:3

        if (!state.fourier_mode)
            @error("fourier_mode required to be true. fourier_mode = $(state.fourier_mode)")
        end

        # Velocity Assignment
        R.IMEX_U[:] .= state.u[:]
        R.IMEX_W[:] .= state.w[:]
        R.IMEX_V[:] .= state.v_F[:]

        # Runge-Kutta Prior Dynamics
        if (i_rk > 1)

            R.IMEX_U[2:ny-1, :, :] +=
                t_stepper.zeta[i_rk] * t_stepper.h_bar[i_rk] * dt * R.RK_U[2:ny-1, :, :]
            R.IMEX_W[2:ny-1, :, :] +=
                t_stepper.zeta[i_rk] * t_stepper.h_bar[i_rk] * dt * R.RK_W[2:ny-1, :, :]

            R.IMEX_V[2:ny_F-1, :, :] +=
                t_stepper.zeta[i_rk] * t_stepper.h_bar[i_rk] * dt * R.RK_V[2:ny_F-1, :, :]

            R.RK_U[:] .= 0.0
            R.RK_V[:] .= 0.0
            R.RK_W[:] .= 0.0

        end

        # Streamwise Pressure-Driven Forcing
        update_flow_attribute!(ns)
        R.IMEX_U[2:ny-1, 1, 1] .+= -dt * t_stepper.h_bar[i_rk] * state.dp0_dx

        # Update Loop
        for iz = 1:domain.nz
            for ix = 1:domain.nx

                # Turbulent Pressure Gradient
                R.IMEX_U[2:ny-1, ix, iz] +=
                    -dt *
                    t_stepper.h_bar[i_rk] *
                    (1.0im * ns.tf.freq_kx[ix]) *
                    state.p[2:ny-1, ix, iz]
                R.IMEX_W[2:ny-1, ix, iz] +=
                    -dt *
                    t_stepper.h_bar[i_rk] *
                    (1.0im * ns.tf.freq_kz[iz]) *
                    state.p[2:ny-1, ix, iz]

                R.IMEX_V[2:ny_F-1, ix, iz] +=
                    -dt * t_stepper.h_bar[i_rk] .*
                    (state.p[2:ny, ix, iz] - state.p[1:ny-1, ix, iz]) ./ domain.dy[:]

                # Homegenous Viscous Dissipation
                freq::ns.FloatType = ns.tf.freq_kx[ix]^2 + ns.tf.freq_kz[iz]^2

                R.RK_U[2:ny-1, ix, iz] = -sim_cond.nu * freq * state.u[2:ny-1, ix, iz]
                R.RK_W[2:ny-1, ix, iz] = -sim_cond.nu * freq * state.w[2:ny-1, ix, iz]

                R.RK_V[2:ny_F-1, ix, iz] = -sim_cond.nu * freq * state.v_F[2:ny_F-1, ix, iz]

            end
        end

        # Transformation
        fourier2physical!(ns.tf, ns.state)

        # Nonlinear Advection - d/dx (UU)
        R.S[:] = state.u[:] .* state.u[:]

        physical2fourier!(ns.tf, R.S, false, true)

        for ix = 1:domain.nx

            R.RK_U[2:ny-1, ix, :] -= 1.0im * ns.tf.freq_kx[ix] * R.S[2:ny-1, ix, :] # d/dx UU

        end

        # Nonlinear Advection - d/dz (WW)
        R.S[:] = state.w[:] .* state.w[:]

        physical2fourier!(ns.tf, R.S, false, true)

        for iz = 1:domain.nz

            R.RK_W[2:ny-1, :, iz] -= 1.0im * ns.tf.freq_kz[iz] * R.S[2:ny-1, :, iz] # d/dz WW

        end

        # Nonlinear Advection - d/dx (UW) & d/dz (UW)
        R.S[:] = state.u[:] .* state.w[:]

        physical2fourier!(ns.tf, R.S, false, true)

        for iz = 1:domain.nz

            R.RK_U[2:ny-1, :, iz] -= 1.0im * ns.tf.freq_kz[iz] * R.S[2:ny-1, :, iz] # d/dz UW

        end

        for ix = 1:domain.nx

            R.RK_W[2:ny-1, ix, :] -= 1.0im * ns.tf.freq_kx[ix] * R.S[2:ny-1, ix, :] # d/dx WU

        end

        # Nonlinear Advection - d/dx (UV)
        for iz = 1:domain.nz
            for ix = 1:domain.nx

                R.S_F[2:ny_F-1, ix, iz] =
                    state.v_F[2:ny_F-1, ix, iz] .* (
                        state.u[2:ny, ix, iz] .* domain.dy_F[2:ny] .+
                        state.u[1:ny-1, ix, iz] .* domain.dy_F[1:ny-1]
                    ) ./ (2.0 * domain.dy[1:ny-1])

            end
        end

        physical2fourier!(ns.tf, R.S_F, true, true)

        for ix = 1:domain.nx

            R.RK_V[2:ny_F-1, ix, :] -= 1.0im * ns.tf.freq_kx[ix] * R.S_F[2:ny_F-1, ix, :]

        end

        # Nonlinear Advection - d/dz (VW)
        for iz = 1:domain.nz
            for ix = 1:domain.nx

                R.S_F[2:ny_F-1, ix, iz] =
                    state.v_F[2:ny_F-1, ix, iz] .* (
                        state.w[2:ny, ix, iz] .* domain.dy_F[2:ny] .+
                        state.w[1:ny-1, ix, iz] .* domain.dy_F[1:ny-1]
                    ) ./ (2.0 * domain.dy[1:ny-1])

            end
        end

        physical2fourier!(ns.tf, R.S_F, true, true)

        for iz = 1:domain.nz

            R.RK_V[2:ny_F-1, :, iz] -= 1.0im * ns.tf.freq_kz[iz] * R.S_F[2:ny_F-1, :, iz]

        end

        # Nonlinear Advection - d/dy (UV)
        for iz = 1:domain.nz
            for ix = 1:domain.nx

                R.S[2:ny-1, ix, iz] =
                    (
                        (state.u[3:ny, ix, iz] + state.u[2:ny-1, ix, iz]) .*
                        state.v_F[3:ny_F-1, ix, iz] -
                        (state.u[2:ny-1, ix, iz] + state.u[1:ny-2, ix, iz]) .*
                        state.v_F[2:ny_F-2, ix, iz]
                    ) ./ (2 * domain.dy_F[2:ny-1])

            end
        end

        physical2fourier!(ns.tf, R.S, false, true)

        R.RK_U[2:ny-1, :, :] -= R.S[2:ny-1, :, :]

        # Nonlinear Advection - d/dy (VW)
        for iz = 1:domain.nz
            for ix = 1:domain.nx

                R.S[2:ny-1, ix, iz] =
                    (
                        (state.w[3:ny, ix, iz] + state.w[2:ny-1, ix, iz]) .*
                        state.v_F[3:ny_F-1, ix, iz] -
                        (state.w[2:ny-1, ix, iz] + state.w[1:ny-2, ix, iz]) .*
                        state.v_F[2:ny_F-2, ix, iz]
                    ) ./ (2 * domain.dy_F[2:ny-1])

            end
        end

        physical2fourier!(ns.tf, R.S, false, true)

        R.RK_W[2:ny-1, :, :] -= R.S[2:ny-1, :, :]

        # Nonlinear Advection - d/dy (VV)
        for iz = 1:domain.nz
            for ix = 1:domain.nx

                R.S_F[2:ny_F-1, ix, iz] =
                    (
                        (state.v_F[3:ny_F, ix, iz] + state.v_F[2:ny_F-1, ix, iz]) .^ 2 -
                        (state.v_F[2:ny_F-1, ix, iz] + state.v_F[1:ny_F-2, ix, iz]) .^ 2
                    ) ./ (4.0 * domain.dy[:])

            end
        end

        physical2fourier!(ns.tf, R.S_F, true, true)

        R.RK_V[2:ny_F-1] -= R.S_F[2:ny_F-1]

        # Assign Runge-Kutta Dynamics to IMEX Dynamics
        R.IMEX_U[2:ny-1] +=
            dt * t_stepper.h_bar[i_rk] * t_stepper.beta[i_rk] * R.RK_U[2:ny-1]
        R.IMEX_W[2:ny-1] +=
            dt * t_stepper.h_bar[i_rk] * t_stepper.beta[i_rk] * R.RK_W[2:ny-1]

        R.IMEX_V[2:ny_F-1] +=
            dt * t_stepper.h_bar[i_rk] * t_stepper.beta[i_rk] * R.RK_V[2:ny_F-1]

        # Transformation Operation for IMEX Dynamics
        fourier2physical!(ns.tf, R.IMEX_U, false)
        fourier2physical!(ns.tf, R.IMEX_W, false)

        fourier2physical!(ns.tf, R.IMEX_V, true)

        for iz = 1:domain.nz
            for ix = 1:domain.nx

                R.IMEX_U[2:ny-1, ix, iz] +=
                    (sim_cond.nu * dt * t_stepper.h_bar[i_rk] / 2.0) * (
                        (
                            (state.u[3:ny, ix, iz] - state.u[2:ny-1, ix, iz]) ./
                            domain.dy[2:ny-1]
                        ) - (
                            (state.u[2:ny-1, ix, iz] - state.u[1:ny-2, ix, iz]) ./
                            domain.dy[1:ny-2]
                        )
                    ) ./ domain.dy_F[2:ny-1]

                R.IMEX_W[2:ny-1, ix, iz] +=
                    (sim_cond.nu * dt * t_stepper.h_bar[i_rk] / 2.0) * (
                        (
                            (state.w[3:ny, ix, iz] - state.w[2:ny-1, ix, iz]) ./
                            domain.dy[2:ny-1]
                        ) - (
                            (state.w[2:ny-1, ix, iz] - state.w[1:ny-2, ix, iz]) ./
                            domain.dy[1:ny-2]
                        )
                    ) ./ domain.dy_F[2:ny-1]


                R.IMEX_V[2:ny_F-1, ix, iz] +=
                    (sim_cond.nu * dt * t_stepper.h_bar[i_rk] / 2.0) * (
                        (
                            (state.v_F[3:ny_F, ix, iz] - state.v_F[2:ny_F-1, ix, iz]) ./
                            domain.dy_F[2:ny]
                        ) ./ domain.dy[:] -
                        (
                            (state.v_F[2:ny_F-1, ix, iz] - state.v_F[1:ny_F-2, ix, iz]) ./
                            domain.dy_F[1:ny-1]
                        ) ./ domain.dy[:]
                    )

            end
        end

        # Solve Intermediate Wall-Normal Velocity
        sol_v::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny_F)

        define_sys_inter_v!(
            ns.solver.sys_inter_v,
            domain,
            sim_cond,
            dt,
            t_stepper.h_bar[i_rk],
        )

        for iz = 1:domain.nz
            for ix = 1:domain.nx

                solve_sys_inter_v!(
                    ns.solver.sys_inter_v,
                    real(R.IMEX_V[:, ix, iz]),
                    sol_v,
                    sim_cond,
                )

                state.v_F[2:ny_F-1, ix, iz] = sol_v[2:ny_F-1]

            end
        end

        # Solve Intermediate Homegeneous Velocity
        sol_u::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)
        sol_w::Vector{ns.FloatType} = zeros(ns.FloatType, domain.ny)

        define_sys_inter_uw!(
            ns.solver.sys_inter_uw,
            domain,
            sim_cond,
            dt,
            t_stepper.h_bar[i_rk],
        )

        for iz = 1:domain.nz
            for ix = 1:domain.nx

                solve_sys_inter_uw!(
                    ns.solver.sys_inter_uw,
                    real(R.IMEX_U[:, ix, iz]),
                    sol_u,
                    sim_cond,
                    'u',
                )
                solve_sys_inter_uw!(
                    ns.solver.sys_inter_uw,
                    real(R.IMEX_W[:, ix, iz]),
                    sol_w,
                    sim_cond,
                    'w',
                )

                state.u[2:ny-1, ix, iz] = sol_u[2:ny-1]
                state.w[2:ny-1, ix, iz] = sol_w[2:ny-1]

            end
        end

        # Transformation
        physical2fourier!(ns.tf, ns.state)

        # Pressure Solve
        continuity_enforcement!(ns, i_rk)

        # println("-- IMEX -- (END) $(i_rk)")
        # println("U | ", maximum(abs.(R.IMEX_U)), " | ", argmax(abs.(R.IMEX_U)))
        # println("V | ", maximum(abs.(R.IMEX_V)), " | ", argmax(abs.(R.IMEX_V)))
        # println("W | ", maximum(abs.(R.IMEX_W)), " | ", argmax(abs.(R.IMEX_W)))
        # println(" ")

        # println("-- Velocity -- (END) $(i_rk)")
        # println("U | ", maximum(abs.(state.u)), " | ", argmax(abs.(state.u)))
        # println("V | ", maximum(abs.(state.v_F)), " | ", argmax(abs.(state.v_F)))
        # println("W | ", maximum(abs.(state.w)), " | ", argmax(abs.(state.w)))
        # println("P | ", maximum(abs.(state.p)), " | ", argmax(abs.(state.p)))
        # println(" ")

    end

    return nothing

end

# Divergence-Free Velocity Pressure Update
function continuity_enforcement!(ns::NavierStokesPropagator, i_rk::Int64)

    if (!ns.state.fourier_mode)

        @error(
            "pressure_update function call with invalid fourier_mode = $(ns.state.fourier_mode)."
        )
        return

    end

    # Simulation Construct
    domain = ns.domain
    state = ns.state

    ny::Int64 = domain.ny
    ny_F::Int64 = domain.ny_F

    # Pressure Solve
    sol_p::Vector{Complex{ns.FloatType}} = zeros(Complex{ns.FloatType}, domain.ny)

    for iz = 1:domain.nz
        for ix = 1:domain.nx

            # Frequency
            kx::ns.FloatType = ns.tf.freq_kx[ix]
            kz::ns.FloatType = ns.tf.freq_kz[iz]

            # Define System
            define_sys_pressure!(ns.solver.sys_pressure, domain, kx, kz)

            # Solve
            solve_sys_pressure_update!(
                ns.solver.sys_pressure,
                sol_p,
                domain,
                kx,
                kz,
                state.u[:, ix, iz],
                state.v_F[:, ix, iz],
                state.w[:, ix, iz],
            )

            # Velocity Update
            state.u[2:ny-1, ix, iz] += -1.0im * kx * sol_p[2:ny-1]
            state.w[2:ny-1, ix, iz] += -1.0im * kz * sol_p[2:ny-1]

            state.v_F[2:ny_F-1, ix, iz] += -(sol_p[2:ny] - sol_p[1:ny-1]) ./ domain.dy[:]

            # Pressure Update
            if (i_rk > 0)
                state.p[:, ix, iz] +=
                    sol_p[:] / (ns.t_stepper.dt * ns.t_stepper.h_bar[i_rk])
            end

        end
    end

    return nothing

end

# Pressure Poisson Solver
function pressure_poisson_solve(ns::NavierStokesPropagator)

    # Simulation Construct
    domain::DomainDescriptor = ns.domain
    state::State = ns.state

    # Domain Grid Size
    nx::Int64 = domain.nx
    nz::Int64 = domain.nz
    ny::Int64 = domain.ny
    ny_F::Int64 = domain.ny_F

    grid_base::NTuple{3,Int64} = (ny, nx, nz)

    # Dynamics Allocation
    R_U::Array{Complex{ns.FloatType},3} = zeros(Complex{ns.FloatType}, grid_base)
    R_V::Array{Complex{ns.FloatType},3} = zeros(Complex{ns.FloatType}, grid_base)
    R_W::Array{Complex{ns.FloatType},3} = zeros(Complex{ns.FloatType}, grid_base)

    S::Array{Complex{ns.FloatType},3} = zeros(Complex{ns.FloatType}, grid_base)

    # Gradient du_i/dx_i
    for iz = 1:domain.nz
        for ix = 1:domain.nx

            R_U[2:ny-1, ix, iz] = 1.0im * ns.tf.freq_kx[ix] * state.u[2:ny-1, ix, iz]
            R_W[2:ny-1, ix, iz] = 1.0im * ns.tf.freq_kz[iz] * state.w[2:ny-1, ix, iz]

            R_V[2:ny-1, ix, iz] =
                (state.v_F[3:ny_F-1, ix, iz] - state.v_F[2:ny_F-2, ix, iz]) ./
                domain.dy_F[2:ny-1]

        end
    end

    fourier2physical!(ns.tf, R_U, false)
    fourier2physical!(ns.tf, R_W, false)
    fourier2physical!(ns.tf, R_V, false)

    # Gradient (du_i/dx_i) (du_i/dx_i)
    R_U[:] .*= R_U[:]
    R_W[:] .*= R_W[:]

    R_V[:] .*= R_V[:]

    physical2fourier!(ns.tf, R_U, false, false)
    physical2fourier!(ns.tf, R_W, false, false)
    physical2fourier!(ns.tf, R_V, false, false)

    # Diagonal Terms
    S[:] = R_U[:] + R_V[:] + R_W[:]

    # Off-Diagonal (du/dy) (dv/dx)
    for iz = 1:domain.nz
        for ix = 1:domain.nx

            R_U[2:ny-1, ix, iz] =
                (state.u[3:ny, ix, iz] - state.u[1:ny-2, ix, iz]) ./
                (2 * domain.dy_F[2:ny-1])

            R_V[2:ny-1, ix, iz] =
                1.0im *
                ns.tf.freq_kx[ix] *
                (0.5 * (state.v_F[3:ny_F-1, ix, iz] + state.v_F[2:ny_F-2, ix, iz]))

        end
    end

    R_V[1, :, :] .= 0.0
    R_V[end, :, :] .= 0.0

    fourier2physical!(ns.tf, R_U, false)
    fourier2physical!(ns.tf, R_V, false)

    R_U[:] = 2.0 * R_U[:] .* R_V[:]

    physical2fourier!(ns.tf, R_U, false, false)

    S[:] .+= R_U[:]

    # Off-Diagonal (dw/dy) (dv/dz)
    for iz = 1:domain.nz
        for ix = 1:domain.nx

            R_W[2:ny-1, ix, iz] =
                (state.w[3:ny, ix, iz] - state.w[1:ny-2, ix, iz]) ./
                (2 * domain.dy_F[2:ny-1])

            R_V[2:ny-1, ix, iz] =
                1.0im *
                ns.tf.freq_kz[iz] *
                (0.5 * (state.v_F[3:ny_F-1, ix, iz] + state.v_F[2:ny_F-2, ix, iz]))

        end
    end

    fourier2physical!(ns.tf, R_W, false)
    fourier2physical!(ns.tf, R_V, false)

    R_W[:] = 2.0 * R_V[:] .* R_W[:]

    physical2fourier!(ns.tf, R_W, false, false)

    S[:] .+= R_W[:]

    # Off-Diagonal (du/dz) (dw/dx)
    for iz = 1:domain.nz
        for ix = 1:domain.nx

            R_U[2:ny-1, ix, iz] = 1.0im * ns.tf.freq_kz[iz] * state.u[2:ny-1, ix, iz]
            R_W[2:ny-1, ix, iz] = 1.0im * ns.tf.freq_kx[ix] * state.w[2:ny-1, ix, iz]

        end
    end

    fourier2physical!(ns.tf, R_U, false)
    fourier2physical!(ns.tf, R_W, false)

    R_U[:] = 2.0 * R_U[:] .* R_W[:]

    physical2fourier!(ns.tf, R_U, false, false)

    S[:] .+= R_U[:]

    # Solve Pressure Poisson Equation
    sol::Vector{Complex{ns.FloatType}} = zeros(Complex{ns.FloatType}, (ny))

    for iz = 1:domain.nz
        for ix = 1:domain.nx

            # Frequency
            kx::ns.FloatType = ns.tf.freq_kx[ix]
            kz::ns.FloatType = ns.tf.freq_kz[iz]

            # Define System
            define_sys_pressure!(ns.solver.sys_pressure, domain, kx, kz)

            # Solve
            solve_sys_pressure_poisson!(ns.solver.sys_pressure, sol, -S[:, ix, iz], kx, kz)

            # Assign
            state.p[:, ix, iz] = sol[:]

        end
    end

end
