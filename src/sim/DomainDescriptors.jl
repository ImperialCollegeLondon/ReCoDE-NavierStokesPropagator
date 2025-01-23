# DomainDescriptors

module DomainDescriptors

#=
Packages
=#

#=
Public API
=#

export DomainDescriptor

#=
Domain Descriptor
=#

struct DomainDescriptor{T<:AbstractFloat}

    # Mesh Size
    nx::Int64
    ny::Int64
    ny_F::Int64
    nz::Int64

    # Domain Length
    Lx::T
    Ly::T
    Lz::T

    # Grid Step
    dx::T
    dz::T

    # Grid Definition
    x::Vector{T}
    z::Vector{T}
    y::Vector{T}
    y_F::Vector{T}

    # Finite Difference Grid
    dy::Vector{T}
    dy_F::Vector{T}

    function DomainDescriptor{T}(
        nx::Int64,
        ny::Int64,
        nz::Int64,
        Lx::T,
        Ly::T,
        Lz::T,
    ) where {T<:AbstractFloat}
        function define_grid!(domain::DomainDescriptor)

            # Periodic Direction
            domain.x[:] =
                range(; start=0.0, stop=domain.Lx, step=domain.dx)[1:(end-1)]
            domain.z[:] =
                range(; start=0.0, stop=domain.Lz, step=domain.dz)[1:(end-1)]

            # Wall-Normal Direction (Base Grid)
            for i in eachindex(domain.y)
                domain.y[domain.ny-(i-1)] = cospi((i - 1) / (domain.ny - 1))
            end

            # Wall-Normal Direction (Fractional Grid)
            for i = 2:(domain.ny_F-1)
                domain.y_F[i] = 0.5 * (domain.y[i] + domain.y[i-1])
            end

            domain.y_F[1] = 2 * domain.y[1] - domain.y_F[2]
            domain.y_F[end] = 2 * domain.y[end] - domain.y_F[end-1]

            # Machine Precision Centerline
            if (mod(domain.ny, 2) != 0)
                domain.y[Int64((domain.ny + 1) / 2)] = 0.0
            end

            if (mod(domain.ny_F, 2) != 0)
                domain.y_F[Int64((domain.ny_F + 1) / 2)] = 0.0
            end

            # Finite Difference Scaling
            for i = 1:(domain.ny-1)
                domain.dy[i] = (domain.y[i+1] - domain.y[i])
            end

            for i = 1:(domain.ny_F-1)
                domain.dy_F[i] = (domain.y_F[i+1] - domain.y_F[i])
            end

            @info("Wall-normal grid definition completed.")

            return nothing
        end

        @info("--- Initialising DomainDescriptor ---")

        # Grid
        ny_F::Int64 = ny + 1

        x::Vector{T} = Vector{T}(undef, nx)
        z::Vector{T} = Vector{T}(undef, nz)
        y::Vector{T} = Vector{T}(undef, ny)
        y_F::Vector{T} = Vector{T}(undef, ny_F)
        dy::Vector{T} = Vector{T}(undef, ny - 1)
        dy_F::Vector{T} = Vector{T}(undef, ny_F - 1)

        # Grid Size
        dx::T = Lx / nx
        dz::T = Lz / nz

        # Initialise
        domain::DomainDescriptor =
            new{T}(nx, ny, ny_F, nz, Lx, Ly, Lz, dx, dz, x, z, y, y_F, dy, dy_F)

        # Define Grid
        define_grid!(domain)

        @info("--- DomainDescriptor initialised ---")
        @info(" ")

        return domain
    end
end

end
