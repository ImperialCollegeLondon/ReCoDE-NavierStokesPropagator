module LinearSolvers

#=
Packages
=#

using ..DomainDescriptors: DomainDescriptor
using ..SimulationConditions: SimulationCondition
using ..States: State

#=
Public API
=#

export TriDiagonalMatrix, LinearSolver
export solve_system_R!, solve_system_C!
export define_sys_inter_v!, solve_sys_inter_v!
export define_sys_inter_uw!, solve_sys_inter_uw!
export define_sys_pressure!, solve_sys_pressure_update!, solve_sys_pressure_poisson!

#=
LinearSolver
=#
"""
Describes a `n x n` matrix which is non-zero in the main, upper and lower off-diagonals.

The first lower diagonal element (TriDiagonalMatrix.L) and the last upper diagonal element (TriDiagonalMatrix.U)
are redundant.
"""
struct TriDiagonalMatrix{T<:Number}

    # Tri-Diagonal Matrix
    L::Vector{T} # Lower Diagonal
    D::Vector{T} # Diagonal
    U::Vector{T} # Upper Diagonal

    function TriDiagonalMatrix{T}(n::Int64) where {T<:Number}

        return new{T}(Vector{T}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n))

    end

end

"""
The LinearSolver encapsulates the required TriDiagonalMatrix system for solving for the velocity and pressure
    update states.
"""
struct LinearSolver{T<:AbstractFloat}

    # Tridiagonal System
    sys_inter_v::TriDiagonalMatrix{T} # Intermediate Wall-Normal Velocity
    sys_inter_uw::TriDiagonalMatrix{T} # Intermediate Homegeneous Velocity
    sys_pressure::TriDiagonalMatrix{Complex{T}} # Pressure Correction

    function LinearSolver{T}(ny::Int64, ny_F::Int64) where {T<:AbstractFloat}

        @info("--- Initialising LinearSolver ---")

        @info("--- LinearSolver initialised ---")
        @info(" ")

        return new(false, TriDiagonalMatrix{T}(ny_F), TriDiagonalMatrix{T}(ny), TriDiagonalMatrix{Complex{T}}(ny))

    end

end

#=
Thomas Algorithm
=#

function solve_system_R!(sys::TriDiagonalMatrix{T}, R::Vector{T}, sol::Vector{T}) where {T<:AbstractFloat}

    # System Dimension
    n::Int64 = length(R)

    # Placeholder Vector
    beta::Vector{T} = Vector{T}(undef, n)
    gamma::Vector{T} = Vector{T}(undef, n)

    # Forward Pass
    beta[1] = sys.D[1]
    gamma[1] = R[1] / beta[1]

    for i = 2:n

        beta[i] = sys.D[i] - sys.L[i] * (sys.U[i-1] / beta[i-1])
        gamma[i] = (-sys.L[i] * gamma[i-1] + R[i]) / beta[i]

    end

    # Backward Pass
    sol[end] = gamma[end]

    for i = n-1:-1:1

        sol[i] = gamma[i] - sol[i+1] * (sys.U[i] / beta[i])

    end

    return nothing

end

function solve_system_C!(sys::TriDiagonalMatrix{Complex{T}}, R::Vector{Complex{T}}, sol::Vector{Complex{T}}) where {T<:AbstractFloat}

    # System Dimension
    n::Int64 = length(R)

    # Placeholder Vector
    beta::Vector{Complex{T}} = Vector{Complex{T}}(undef, n)
    gamma::Vector{Complex{T}} = Vector{Complex{T}}(undef, n)

    # Forward Pass
    beta[1] = sys.D[1]
    gamma[1] = R[1] / beta[1]

    for i = 2:n

        beta[i] = sys.D[i] - sys.L[i] * (sys.U[i-1] / beta[i-1])
        gamma[i] = (-sys.L[i] * gamma[i-1] + R[i]) / beta[i]

    end

    # Backward Pass
    sol[end] = gamma[end]

    for i = n-1:-1:1

        sol[i] = gamma[i] - sol[i+1] * (sys.U[i] / beta[i])

    end

    return nothing

end

#=
Intermediate Runge-Kutta Wall-Normal Velocity
=#

# System Definition
function define_sys_inter_v!(sys::TriDiagonalMatrix{T}, domain::DomainDescriptor, sim_cond::SimulationCondition, dt::T, h_bar::T) where {T<:AbstractFloat}

    # Scaling Factor
    scaling::T = 0.5 * sim_cond.nu * dt * h_bar

    # Lower Diagonal
    sys.L[2:end-1] = -scaling ./ (domain.dy_F[1:end-1] .* domain.dy[:])

    # Main Diagonal
    sys.D[2:end-1] = 1.0 .+ scaling ./ (domain.dy_F[2:end] .* domain.dy[:]) + scaling ./ (domain.dy_F[1:end-1] .* domain.dy[:])

    # Upper Diagonal
    sys.U[2:end-1] = -scaling ./ (domain.dy_F[2:end] .* domain.dy[:])

end

# Solve System
function solve_sys_inter_v!(sys::TriDiagonalMatrix{T}, R::Vector{T}, sol::Vector{T}, sim_cond::SimulationCondition) where {T<:AbstractFloat}

    # Apply Boundary Condition (Lower Wall)
    sys.D[1] = 1.0
    sys.L[1] = 0.0
    sys.U[1] = 0.0
    R[1] = sim_cond.bound_cond[(vel='v', wall='l')]

    # Apply Boundary Condition (Upper Wall)
    sys.D[end] = 1.0
    sys.L[end] = 0.0
    sys.U[end] = 0.0
    R[end] = sim_cond.bound_cond[(vel='v', wall='u')]

    # Solve System
    solve_system_R!(sys, R, sol)

end

#=
Intermediate Runge-Kutta Homegeneous Velocity
=#

# System Definition
function define_sys_inter_uw!(sys::TriDiagonalMatrix{T}, domain::DomainDescriptor, sim_cond::SimulationCondition, dt::T, h_bar::T) where {T<:AbstractFloat}

    # Scaling Factor
    scaling::T = 0.5 * sim_cond.nu * dt * h_bar

    # Lower Diagonal
    sys.L[2:end-1] = -scaling ./ (domain.dy_F[2:end-1] .* domain.dy[1:end-1])

    # Main Diagonal
    sys.D[2:end-1] = 1.0 .+ scaling ./ (domain.dy_F[2:end-1] .* domain.dy[2:end]) + scaling ./ (domain.dy_F[2:end-1] .* domain.dy[1:end-1])

    # Upper Diagonal
    sys.U[2:end-1] = -scaling ./ (domain.dy_F[2:end-1] .* domain.dy[2:end])

end

# Solve System
function solve_sys_inter_uw!(sys::TriDiagonalMatrix{T}, R::Vector{T}, sol::Vector{T}, sim_cond::SimulationCondition, vel_comp::Char) where {T<:AbstractFloat}

    # Apply Boundary Condition (Lower Wall)
    sys.D[1] = 1.0
    sys.L[1] = 0.0
    sys.U[1] = 0.0
    R[1] = sim_cond.bound_cond[(vel=vel_comp, wall='l')]

    # Apply Boundary Condition (Upper Wall)
    sys.D[end] = 1.0
    sys.L[end] = 0.0
    sys.U[end] = 0.0
    R[end] = sim_cond.bound_cond[(vel=vel_comp, wall='u')]

    # Solve System
    solve_system_R!(sys, R, sol)

end

#=
Pressure Correction Poisson Solve
=#

# System Definition
function define_sys_pressure!(sys::TriDiagonalMatrix{Complex{T}}, domain::DomainDescriptor, kx::T, kz::T) where {T<:AbstractFloat}

    # Lower Diagonal
    sys.L[2:end-1] = 1.0 ./ (domain.dy_F[2:end-1] .* domain.dy[1:end-1])

    # Main Diagonal
    sys.D[2:end-1] = -(kx^2 + kz^2) .- (1.0 ./ (domain.dy_F[2:end-1] .* domain.dy[1:end-1])) - (1.0 ./ (domain.dy_F[2:end-1] .* domain.dy[2:end]))

    # Upper Diagonal
    sys.U[2:end-1] = 1.0 ./ (domain.dy_F[2:end-1] .* domain.dy[2:end])

end

# Solve System
function solve_sys_pressure_update!(sys::TriDiagonalMatrix{Complex{T}}, sol::Vector{Complex{T}}, domain::DomainDescriptor,
    kx::T, kz::T, u::Vector{Complex{T}}, v_F::Vector{Complex{T}}, w::Vector{Complex{T}}) where {T<:AbstractFloat}

    # Dynamics Vector
    R::Vector{Complex{T}} = Vector{Complex{T}}(undef, domain.ny)

    R[:] = 1.0im * kx * u[:] + 1.0im * kz * w[:] + (v_F[2:end] - v_F[1:end-1]) ./ domain.dy_F[:]

    # Boundary Condition
    if (kx == 0) && (kz == 0)

        # Apply Boundary Condition (Lower Wall)
        sys.D[1] = 1.0
        sys.L[1] = 0.0
        sys.U[1] = 0.0
        R[1] = 0.0

    else

        # Apply Boundary Condition (Lower Wall)
        sys.D[1] = 1.0
        sys.L[1] = 0.0
        sys.U[1] = -1.0
        R[1] = 0.0

    end

    # Apply Boundary Condition (Upper Wall)
    sys.D[end] = -1.0
    sys.L[end] = 1.0
    sys.U[end] = 0.0
    R[end] = 0.0

    # Solve System
    solve_system_C!(sys, R, sol)

end

# Solve System
function solve_sys_pressure_poisson!(sys::TriDiagonalMatrix{Complex{T}}, sol::Vector{Complex{T}}, R::Vector{Complex{T}},
    kx::T, kz::T,) where {T<:AbstractFloat}

    # Boundary Condition
    if (kx == 0.0) && (kz == 0.0)

        # Apply Boundary Condition (Lower Wall)
        sys.D[1] = 1.0
        sys.L[1] = 0.0
        sys.U[1] = 0.0
        R[1] = 0.0

    else

        # Apply Boundary Condition (Lower Wall)
        sys.D[1] = 1.0
        sys.L[1] = 0.0
        sys.U[1] = -1.0
        R[1] = 0.0

    end

    # Apply Boundary Condition (Upper Wall)
    sys.D[end] = -1.0
    sys.L[end] = 1.0
    sys.U[end] = 0.0
    R[end] = 0.0

    # Solve System
    solve_system_C!(sys, R, sol)

end

end
