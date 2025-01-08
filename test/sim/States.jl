# Test for State

#=
Test Set
=#

@testset "State Array Size" begin

    # Initialise state
    state::State = State{Float64}(10.0, 10, 10, 10)

    # Array Size Test
    @test size(state.u) == (10, 10, 10)
    @test size(state.v) == (10, 10, 10)
    @test size(state.w) == (10, 10, 10)
    @test size(state.v_F) == (11, 10, 10)

    @test length(state.u_bar) == 10
    @test length(state.v_bar) == 10
    @test length(state.w_bar) == 10
    @test length(state.uu_bar) == 10
    @test length(state.uv_bar) == 10
    @test length(state.uw_bar) == 10
    @test length(state.vv_bar) == 10
    @test length(state.vw_bar) == 10
    @test length(state.ww_bar) == 10
end
