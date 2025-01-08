# Test for Transformation

#=
Test Set
=#

@testset "Transformation Struct Property" begin

    # Initialise Simulation Construct
    Lx::Float64, Lz::Float64 = (10.0, 10.0)
    nx::Int64, ny::Int64, nz::Int64 = (12, 10, 12)

    domain::DomainDescriptor = DomainDescriptor{Float64}(nx, ny, nz, Lx, 2.0, Lz)
    state::State = State{Float64}(10.0, nx, ny, nz)

    # Initialise Transformation
    tf::Transformer = Transformer{Float64}(
        "test/fftw_wisdom", state, domain; wisdom_flag="exhaustive", test_mode=true
    )

    @test length(tf.freq_kx) == nx
    @test length(tf.freq_kz) == nz

    @test tf.freq_kx_max == 4 * 2pi / Lx
    @test tf.freq_kz_max == 4 * 2pi / Lz
end