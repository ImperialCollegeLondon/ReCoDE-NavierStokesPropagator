# State

module States

#=
Packages
=#

#=
Public API
=#

export State

#=
State
=#

mutable struct State{T<:AbstractFloat}

    # Simulation Time
    t_sim::T # Current Time
    T_sim::T # Total Time
    nt::Int64 # Current Time Step

    # Velocity-Pressure State
    u::Array{Complex{T},3}
    v::Array{Complex{T},3}
    v_F::Array{Complex{T},3}
    w::Array{Complex{T},3}
    p::Array{Complex{T},3}
    fourier_mode::Bool
    frac_mode::Bool

    # Statistical State
    n_sample::Int64
    u_bar::Vector{T}
    v_bar::Vector{T}
    w_bar::Vector{T}
    uu_bar::Vector{T}
    uv_bar::Vector{T}
    uw_bar::Vector{T}
    vv_bar::Vector{T}
    vw_bar::Vector{T}
    ww_bar::Vector{T}

    # Flow Attribute
    dp0_dx::T # Mean Pressure Gradient
    u_bulk::T # Bulk Velocity
    u_tau_upr::T # Upper Wall Shear Stress
    u_tau_lwr::T # Lower Wall Shear Stress

    function State{T}(T_sim::T, nx::Int64, ny::Int64, nz::Int64) where {T<:AbstractFloat}
        @info("--- Initialising State ---")
        @info("--- State initialised ---")
        @info(" ")

        return new{T}(
            0.0,
            T_sim,
            0,
            zeros(Complex{T}, (ny, nx, nz)), # Velocity-Pressure State
            zeros(Complex{T}, (ny, nx, nz)),
            zeros(Complex{T}, (ny + 1, nx, nz)),
            zeros(Complex{T}, (ny, nx, nz)),
            zeros(Complex{T}, (ny, nx, nz)),
            false,
            false,
            0, # Statistical State
            zeros(T, (ny)),
            zeros(T, (ny)),
            zeros(T, (ny)),
            zeros(T, (ny)),
            zeros(T, (ny)),
            zeros(T, (ny)),
            zeros(T, (ny)),
            zeros(T, (ny)),
            zeros(T, (ny)),
            0.0,
            0.0,
            0.0,
            0.0,
        ) # Flow Attribute
    end
end

end
