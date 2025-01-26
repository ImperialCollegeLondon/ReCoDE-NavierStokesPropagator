# Transformations

module Transformations

#=
Package
=#

using FFTW

using ..States: State
using ..DomainDescriptors: DomainDescriptor

#=
Public API
=#

export Transformer
export y2y_F!, y_F2y!
export physical2fourier!, fourier2physical!

#=
Transformer
=#

# FFTW Wisdom Flag
wisdom_dict::Dict{String,UInt32} = Dict()

wisdom_dict["estimate"] = FFTW.ESTIMATE
wisdom_dict["measure"] = FFTW.MEASURE
wisdom_dict["patient"] = FFTW.PATIENT
wisdom_dict["exhaustive"] = FFTW.EXHAUSTIVE

struct Transformer{T<:AbstractFloat}

    # FFTW Plan (Base Grid)
    fft_plan::Any
    ifft_plan::Any

    # FFTW Plan (Fractional Grid)
    fft_plan_F::Any
    ifft_plan_F::Any

    # FFTW Frequency
    freq_kx::Vector{T}
    freq_kz::Vector{T}

    # De-aliasing Frequency Bound
    freq_kx_max::T
    freq_kz_max::T

    function Transformer{T}(
        dir_wisdom::String,
        state::State{T},
        domain::DomainDescriptor;
        wisdom_flag::String="measure",
        reset_wisdom::Bool=false,
        test_mode::Bool=false,
    ) where {T<:AbstractFloat}
        @info("--- Initialising Transformer ---")

        # Directory Check
        if (!test_mode)
            if (!isdir(dir_wisdom))
                mkdir(dir_wisdom)
            end

            if (isfile(string(dir_wisdom, "/fft"))) && (!reset_wisdom)
                @info("Loading FFTW wisdom from $(dir_wisdom).")
                FFTW.import_wisdom(string(dir_wisdom, "/fft"))
            end
        end

        @info("Initiating FFTW wisdom planning.")

        # Planning
        @time begin
            fft_plan = plan_fft!(state.u, (2, 3); flags=wisdom_dict[wisdom_flag])
            ifft_plan = plan_bfft!(state.u, (2, 3); flags=wisdom_dict[wisdom_flag])

            fft_plan_F = plan_fft!(state.v_F, (2, 3); flags=wisdom_dict[wisdom_flag])
            ifft_plan_F = plan_bfft!(state.v_F, (2, 3); flags=wisdom_dict[wisdom_flag])
        end

        if (!test_mode)
            FFTW.export_wisdom(string(dir_wisdom, "/fft"))
        end

        @info("FFTW wisdom planning complete.")

        # Frequency
        n::Int64 = size(state.u, 2)
        freq_kx::Vector{T} = fftfreq(n) * n * (2pi / domain.Lx)
        freq_kx_max::T = floor(n / 3) * (2pi / domain.Lx)

        n = size(state.u, 3)
        freq_kz::Vector{T} = fftfreq(n) * n * (2pi / domain.Lz)
        freq_kz_max::T = floor(n / 3) * (2pi / domain.Lz)

        @info("--- Transformer initialised ---")
        @info(" ")

        return new{T}(
            fft_plan,
            ifft_plan,
            fft_plan_F,
            ifft_plan_F,
            freq_kx,
            freq_kz,
            freq_kx_max,
            freq_kz_max,
        )
    end
end

#=
FFT Transformation
=#

"""
Converts the physical representation of State arrays to frequency domain via FFT.
"""
function physical2fourier!(tf::Transformer{T}, state::State{T}) where {T<:AbstractFloat}

    # Mode Check
    if (state.fourier_mode)
        @debug(
            "Applying a physical2fourier! fftw transformation when fourier_mode = $(state.fourier_mode)."
        )
        return nothing
    end

    # FFT
    tf.fft_plan * state.u
    tf.fft_plan * state.v
    tf.fft_plan * state.w
    tf.fft_plan * state.p

    tf.fft_plan_F * state.v_F

    # Normalisation
    n::T = prod(size(state.u)[2:3])

    state.u[:] ./= n
    state.v[:] ./= n
    state.w[:] ./= n
    state.p[:] ./= n

    state.v_F[:] ./= n

    # Update Mode
    state.fourier_mode = true

    return nothing
end

"""
Converts the physical representation of State arrays to frequency domain via FFT for a specfic array.
"""
function physical2fourier!(
    tf::Transformer{T},
    arr::Array{Complex{T},3},
    frac_mode::Bool,
    dealias_mode::Bool,
) where {T<:AbstractFloat}

    # FFT
    if (frac_mode)
        tf.fft_plan_F * arr

    else
        tf.fft_plan * arr
    end

    # De-aliasing
    if (dealias_mode)
        for ix = 1:nx
            if (abs(tf.freq_kx[ix]) > tf.freq_kx_max)
                arr[:, ix, :] .= 0.0 + 0.0im
            end
        end

        for iz = 1:nz
            if (abs(tf.freq_kz[iz]) > tf.freq_kz_max)
                arr[:, :, iz] .= 0.0 + 0.0im
            end
        end
    end

    # Normalisation
    n::T = prod(size(arr)[2:3])

    arr[:] ./= n

    return nothing
end

"""
Converts the frequency domain representation of State arrays to physical domain via iFFT.
"""
function fourier2physical!(tf::Transformer{T}, state::State{T}) where {T<:AbstractFloat}

    # Mode Check
    if (!state.fourier_mode)
        @debug(
            "Applying a fourier2physical! fftw transformation when fourier_mode = $(state.fourier_mode)."
        )
        return nothing
    end

    # IFFT
    tf.ifft_plan * state.u
    tf.ifft_plan * state.v
    tf.ifft_plan * state.w
    tf.ifft_plan * state.p

    tf.ifft_plan_F * state.v_F

    # Update Mode
    state.fourier_mode = false

    state.u[:] = real(state.u[:])
    state.v[:] = real(state.v[:])
    state.w[:] = real(state.w[:])
    state.p[:] = real(state.p[:])

    state.v_F[:] = real(state.v_F[:])

    return nothing
end

"""
Converts the frequency domain representation of State arrays to physical domain via iFFT for a specfic array.
"""
function fourier2physical!(
    tf::Transformer{T},
    arr::Array{Complex{T},3},
    frac_mode::Bool,
) where {T<:AbstractFloat}

    # IFFT
    if (frac_mode)
        tf.ifft_plan_F * arr

    else
        tf.ifft_plan * arr
    end

    arr[:] = real(arr[:])

    return nothing
end

#=
Grid Transformation
=#

"""
Transforms the grid representation of wall-normal velocity from physical to fractional grid.
"""
function y2y_F!(state::State{T}) where {T<:AbstractFloat}
    if (state.frac_mode)
        @debug("Applying a y2y_F! grid transformation when frac_mode = $(state.frac_mode).")
        return nothing
    end

    # Array
    arr = state.v
    arr_F = state.v_F

    # Sizing
    ny_F::Int64 = size(arr_F, 1)

    # Transformation
    arr_F[1, :, :] .= arr[1, :, :]
    arr_F[end, :, :] .= arr[end, :, :]

    for iy = 2:(ny_F-1)
        arr_F[iy, :, :] .= 0.5 * (arr[iy, :, :] + arr[iy-1, :, :])
    end

    # Swap Mode
    state.frac_mode = true

    return nothing
end

"""
Transforms the grid representation of wall-normal velocity from fractional to physical grid.
"""
function y_F2y!(state::State)
    if (!state.frac_mode)
        @debug("Applying a y_F2y! grid transformation when frac_mode = $(state.frac_mode).")
        return nothing
    end

    # Array
    arr = state.v
    arr_F = state.v_F

    # Sizing
    n::NTuple{3,Int64} = size(arr)

    ny::Int64 = n[1]

    # Enforce Boundary Condition
    arr[1, :, :] .= arr_F[1, :, :]
    arr[end, :, :] .= arr_F[end, :, :]

    # Transformation Loop
    for iy = 2:(ny-1)
        arr[iy, :, :] = 2 * arr_F[iy, :, :] - arr[iy-1, :, :]
    end

    # Switch Mode
    state.frac_mode = false

    return nothing
end

end
