# Test for DomainDescriptor

#=
Test Set
=#

@testset "DomainDescriptor Array Size" begin

    # Initialise Domain
    domain::DomainDescriptor = DomainDescriptor{Float64}(10, 10, 10, 10.0, 2.0, 10.0)

    # Array Size Test
    @test length(domain.x) == 10
    @test length(domain.z) == 10

    @test length(domain.y) == 10
    @test length(domain.dy) == length(domain.y) - 1

    @test length(domain.y_F) == length(domain.y) + 1
    @test length(domain.dy_F) == length(domain.y_F) - 1
end

@testset "DomainDescriptor Values" begin

    # Initialise Domain
    domain::DomainDescriptor = DomainDescriptor{Float64}(10, 10, 10, 10.0, 2.0, 10.0)

    # Homogeneous Direction
    @test domain.x[1] == 0.0
    @test domain.x[end] + domain.dx == domain.Lx

    @test domain.z[1] == 0.0
    @test domain.z[end] + domain.dz == domain.Lz

    # Wall-Normal Direction
    @test domain.y[1] == -1.0
    @test domain.y[end] == 1.0

    @test domain.y_F[1] < domain.y[1]
    @test domain.y_F[end] > domain.y[end]

    # Centerline Machine Precision
    @test domain.y_F[6] == 0.0

    domain_odd::DomainDescriptor = DomainDescriptor{Float64}(10, 11, 10, 10.0, 2.0, 10.0)
    @test domain_odd.y[6] == 0.0
end