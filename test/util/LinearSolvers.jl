# Test for LinearSolvers

#=
Test Set
=#

@testset "LinearSolvers Test" begin

    # No. of Node
    n::Int64 = 5

    # System of Equation
    sys_R::TriDiagonalMatrix{Float64} = TriDiagonalMatrix{Float64}(n)
    sys_R.L[:] = -1.0 * ones(5)
    sys_R.D[:] = 2.0 * ones(5)
    sys_R.U[:] = -1.0 * ones(5)

    sys_C::TriDiagonalMatrix{ComplexF64} = TriDiagonalMatrix{ComplexF64}(n)
    sys_C.L[:] = -1.0 * ones(5)
    sys_C.D[:] = 2.0 * ones(5)
    sys_C.U[:] = -1.0 * ones(5)

    # Solution
    A::Matrix{Float64} = [2 -1 0 0 0; -1 2 -1 0 0; 0 -1 2 -1 0; 0 0 -1 2 -1; 0 0 0 -1 2]
    R::Vector{Float64} = ones(Float64, (n))
    R_C::Vector{ComplexF64} = ones(ComplexF64, (n))

    x::Vector{Float64} = inv(A) * R

    # Solve Real System
    x_R::Vector{Float64} = Vector{Float64}(undef, (n))
    solve_system_R!(sys_R, R, x_R)

    @test x_R ≈ x

    # Solve Complex System
    x_C::Vector{ComplexF64} = Vector{ComplexF64}(undef, (n))
    solve_system_C!(sys_C, R_C, x_C)

    @test real(x_C) ≈ x

end
