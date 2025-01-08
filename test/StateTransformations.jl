# Test for Transformation

#=
Test Set
=#

@testset "Transformation Public API" begin

    # Initialise Simulation Construct
    Lx::Float64, Lz::Float64 = (10.0, 10.0)
    nx::Int64, ny::Int64, nz::Int64 = (12, 10, 12)

    grid_physical = (ny, nx, nz)
    grid_frac = (ny + 1, nx, nz)

    domain::DomainDescriptor = DomainDescriptor{Float64}(nx, ny, nz, Lx, 2.0, Lz)
    state::State = State{Float64}(10.0, nx, ny, nz)

    # Initialise Transformation
    tf::Transformer = Transformer{Float64}(
        "test/fftw_wisdom",
        state,
        domain;
        wisdom_flag = "exhaustive",
        test_mode = true,
    )

    # Random Array
    state.v[:] = rand(Float64, grid_physical)
    state.v_F[:] = rand(Float64, grid_frac)

    # Placeholder
    pl_v::Array{ComplexF64,3} = Array{Float64,3}(undef, (grid_physical))
    pl_v[:] = state.v[:]

    # FFTW - Physical to Fourier
    physical2fourier!(tf, state)
    @test state.fourier_mode == true

    # FFTW - Fourier to Physical
    fourier2physical!(tf, state)
    @test state.fourier_mode == false

    # Fractional Transformation - Physical to Fractional
    y2y_F!(state)
    @test state.frac_mode == true

    # Fractional Transformation - Fractional to Physical
    y_F2y!(state)
    @test state.frac_mode == false

    # Inversion Check
    @test pl_v[:] ≈ state.v[:]
end
