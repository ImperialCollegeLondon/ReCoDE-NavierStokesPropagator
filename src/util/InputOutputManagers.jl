# InputOutputManager

module InputOutputManagers

#=
Packages
=#

using DocStringExtensions

using HDF5

using ..States: State
using ..DomainDescriptors: DomainDescriptor

#=
Public API
=#

export InputOutputManager
export write_grid
export write_flowfield, read_flowfield!
export write_flow_statistics, write_attribute

#=
InputOutputManager
=#
"""
Contains the path directory for reading and writing flowfield and grid points. Also contains 
    the I/O frequency.
"""
struct InputOutputManager

    # Directory
    input_root::String
    output_root::String

    # Output Frequency
    freq_state::Int64
    freq_stats::Int64
    min_step_stats::Int64

    function InputOutputManager(
        input_root::String,
        output_root::String,
        freq_state::Int64,
        freq_stats::Int64,
        min_step_stats::Int64,
    )

        @info("--- Initialising InputOutputManager ---")

        # Reset Output Directory
        if isdir(output_root)

            chmod(output_root, 0o777; recursive = true)
            rm(output_root; recursive = true)

        end

        mkdir(output_root)

        # Flowfield Directory
        mkdir(string(output_root, "/flowfield"))

        # Statistics File
        mkdir(string(output_root, "/statistics"))

        @info("--- InputOutputManager initialised ---")
        @info(" ")

        return new(input_root, output_root, freq_state, freq_stats, min_step_stats)

    end

end

#=
Write Methods
=#
"""
$(SIGNATURES)

Write the wall-normal grid points to file.
"""
function write_grid(
    io_m::InputOutputManager,
    domain::DomainDescriptor{T},
) where {T<:AbstractFloat}

    @info("Writing domain grid to output directory.")

    # Output Path
    output_path::String = string(io_m.output_root, "/grid.h5")

    # Array Size
    (nx::Int64, nz::Int64, ny::Int64) = (domain.nx, domain.nz, domain.ny)

    # Write Statistics and Attribute
    h5open(output_path, "w") do fid

        g::HDF5.Group = open_group(fid, "/")

        # Attribute
        attrs(g)["mesh_size"] = [nx; ny; nz] # Mesh Size

        # Y-Domain
        dset::HDF5.Dataset =
            create_dataset(g, "y", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        write(dset, domain.y)

        # X-Domain
        dset = create_dataset(g, "x", datatype(T), dataspace((nx,)); chunk = ((nx,)))

        write(dset, domain.x)

        # Z-Domain
        dset = create_dataset(g, "z", datatype(T), dataspace((nz,)); chunk = ((nz,)))

        write(dset, domain.z)

    end

    @info("Domain grid written to output directory.")

    return nothing

end

# Write Flowfield
"""
$(SIGNATURES)

Rearranges the flowfield state to (x,y,z) from (y,x,z) and writes to file.
"""
function write_flowfield(
    io_m::InputOutputManager,
    state::State,
    domain::DomainDescriptor{T},
) where {T<:AbstractFloat}

    @info("Writing state t = $(state.t_sim), nt = $(state.nt) to output directory.")

    # Format Output Path
    output_path::String = string(io_m.output_root, "/flowfield/state_nt")

    numeric_vec::Vector{Int64} = digits(state.nt, base = 10, pad = 8)

    for i = 8:-1:1
        output_path = string(output_path, "$(numeric_vec[i])")
    end

    output_path = string(output_path, ".h5")

    # Array Size
    (nx::Int64, nz::Int64, ny::Int64) = (domain.nx, domain.nz, domain.ny)

    # Write State and Attribute
    h5open(output_path, "w") do fid

        g::HDF5.Group = open_group(fid, "/")

        # Attribute
        attrs(g)["t"] = state.t_sim # Simulation Time
        attrs(g)["mesh_size"] = [nx; ny; nz] # Mesh Size

        # Placeholder Array
        pl_arr::Array{T,3} = Array{T,3}(undef, (nx, ny, nz))

        # U Velocity
        dset::HDF5.Dataset = create_dataset(
            g,
            "u",
            datatype(T),
            dataspace((nx, ny, nz));
            chunk = ((nx, 1, nz)),
        )

        pl_arr[:] = permutedims(real(state.u), (2, 1, 3))
        write(dset, pl_arr)

        # V Velocity
        dset = create_dataset(
            g,
            "v",
            datatype(T),
            dataspace((nx, ny, nz));
            chunk = ((nx, 1, nz)),
        )

        pl_arr[:] = permutedims(real(state.v), (2, 1, 3))
        write(dset, pl_arr)

        # W Velocity
        dset = create_dataset(
            g,
            "w",
            datatype(T),
            dataspace((nx, ny, nz));
            chunk = ((nx, 1, nz)),
        )

        pl_arr[:] = permutedims(real(state.w), (2, 1, 3))
        write(dset, pl_arr)

        # P Velocity
        dset = create_dataset(
            g,
            "p",
            datatype(T),
            dataspace((nx, ny, nz));
            chunk = ((nx, 1, nz)),
        )

        pl_arr[:] = permutedims(real(state.p), (2, 1, 3))
        write(dset, pl_arr)

    end

    @info("State t = $(state.t_sim), nt = $(state.nt) written to output directory.")

    return nothing

end

# Write Flow Statistics
"""
$(SIGNATURES)

Rearranges the statistical state to (x,y,z) from (y,x,z) and writes to file.
"""
function write_flow_statistics(
    io_m::InputOutputManager,
    state::State,
    domain::DomainDescriptor{T},
) where {T<:AbstractFloat}

    @info("Writing statistics t = $(state.t_sim), nt = $(state.nt) to output directory.")

    # Format Output Path
    output_path::String = string(io_m.output_root, "/statistics/statistics_nt")

    numeric_vec::Vector{Int64} = digits(state.nt, base = 10, pad = 8)

    for i = 8:-1:1
        output_path = string(output_path, "$(numeric_vec[i])")
    end

    output_path = string(output_path, ".h5")

    # Array Size
    (nx::Int64, nz::Int64, ny::Int64) = (domain.nx, domain.nz, domain.ny)

    # Write Statistics and Attribute
    h5open(output_path, "w") do fid

        g::HDF5.Group = open_group(fid, "/")

        # Attribute
        attrs(g)["t"] = state.t_sim # Simulation Time
        attrs(g)["n_sample"] = state.n_sample
        attrs(g)["mesh_size"] = [nx; ny; nz] # Mesh Size

        # Placeholder Array
        pl_arr::Vector{T} = Vector{T}(undef, ny)

        # Y-Domain
        dset::HDF5.Dataset =
            create_dataset(g, "y", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(domain.y)
        write(dset, pl_arr)

        # U Velocity
        dset = create_dataset(g, "u_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.u_bar)
        write(dset, pl_arr)

        # V Velocity
        dset = create_dataset(g, "v_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.v_bar)
        write(dset, pl_arr)

        # W Velocity
        dset = create_dataset(g, "w_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.w_bar)
        write(dset, pl_arr)

        # UU Velocity
        dset = create_dataset(g, "uu_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.uu_bar)
        write(dset, pl_arr)

        # UV Velocity
        dset = create_dataset(g, "uv_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.uv_bar)
        write(dset, pl_arr)

        # UW Velocity
        dset = create_dataset(g, "uw_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.uw_bar)
        write(dset, pl_arr)

        # VV Velocity
        dset = create_dataset(g, "vv_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.vv_bar)
        write(dset, pl_arr)

        # VW Velocity
        dset = create_dataset(g, "vw_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.vw_bar)
        write(dset, pl_arr)

        # WW Velocity
        dset = create_dataset(g, "ww_bar", datatype(T), dataspace((ny,)); chunk = ((ny,)))

        pl_arr[:] = real(state.ww_bar)
        write(dset, pl_arr)

    end

    @info("Statistics t = $(state.t_sim), nt = $(state.nt) written to output directory.")

    return nothing

end

# Write Attribute
"""
$(SIGNATURES)

Write the simulation attribute on the fly to a logger.
"""
function write_attribute(
    io_m::InputOutputManager,
    state::State{T},
    dt::T,
    nu::T,
) where {T<:AbstractFloat}

    # Output Path
    output_path::String = string(io_m.output_root, "/flowfield/attribute.h5")

    # # Write Flow Attribute
    # h5open(output_path, "cw") do fid

    #     # Output Group
    #     output_group::String = "/nt"
    #     numeric_vec::Vector{Int64} = digits(state.nt, base=10, pad=8)

    #     for i = 8:-1:1
    #         output_group = string(output_group, "$(numeric_vec[i])")
    #     end

    #     g::HDF5.Group = create_group(fid, output_group)

    #     # Time Attribute
    #     attrs(g)["t"] = state.t_sim
    #     attrs(g)["nt"] = state.nt
    #     attrs(g)["dt"] = dt

    #     # Flow Attribute
    #     attrs(g)["dp0_dx"] = state.dp0_dx
    #     attrs(g)["u_bulk"] = state.u_bulk
    #     attrs(g)["u_tau_lwr"] = state.u_tau_lwr
    #     attrs(g)["u_tau_upr"] = state.u_tau_upr
    #     attrs(g)["re_tau_lwr"] = state.u_tau_lwr / nu
    #     attrs(g)["re_tau_upr"] = state.u_tau_upr / nu

    # end

    @info("---")
    @info(" ")

    @info(" t = $(state.t_sim) | nt = $(state.nt) | dt = $(dt)")
    @info(" dp0_dx = $(state.dp0_dx) | u_bulk = $(state.u_bulk)")
    @info(" Lower Wall: u_tau = $(state.u_tau_lwr) | Re_tau = $(state.u_tau_lwr / nu)")
    @info(" Upper Wall: u_tau = $(state.u_tau_upr) | Re_tau = $(state.u_tau_upr / nu)")

    @info(" ")
    @info("---")

    return nothing

end


#=
Read Methods
=#

# Read Flowfield
"""
$(SIGNATURES)

Read the flowfield to a file and convert from (x,y,z) grid to the (y,x,z) grid.
"""
function read_flowfield!(io_m::InputOutputManager, state::State)

    # Format Input Path
    input_path::String = string(io_m.input_root, "/start.h5")

    @info("Reading initial condition ($(input_path)).")

    # Mesh Size Compatiblity Check
    h5open(input_path, "r") do fid

        g::HDF5.Group = open_group(fid, "/")

        (nx::Int64, ny::Int64, nz::Int64) = read_attribute(g, "mesh_size")

        # Array Size
        (ny_arr::Int64, nx_arr::Int64, nz_arr::Int64) = size(state.u)

        # Compatiblity Check
        if ((nx != nx_arr) | (nz != nz_arr) | (ny != ny_arr))

            throw(
                DimensionMismatch(
                    "Mismatching mesh size of start.h5 ($(input_path)) with simulation mesh.",
                ),
            )

        end

    end

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

    @info("Initial condition read completed.")

    return nothing

end

end
