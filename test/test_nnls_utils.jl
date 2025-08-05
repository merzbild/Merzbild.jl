@testset "nnls_merging utils" begin
    base_mim = Merzbild.base_multi_index_moments()
    ref_mom_indices = [[0, 0, 0],
                       [1, 0, 0], [0, 1, 0], [0, 0, 1],
                       [2, 0, 0], [0, 2, 0], [0, 0, 2]]

    @test length(base_mim) == length(ref_mom_indices)
    for ref_index in ref_mom_indices
        @test (ref_index in base_mim) == true
    end

    @test Merzbild.vx_sign(1) == -1
    @test Merzbild.vy_sign(1) == -1
    @test Merzbild.vz_sign(1) == -1
    
    @test Merzbild.vx_sign(2) == 1
    @test Merzbild.vy_sign(2) == -1
    @test Merzbild.vz_sign(2) == -1

    @test Merzbild.vx_sign(3) == -1
    @test Merzbild.vy_sign(3) == 1
    @test Merzbild.vz_sign(3) == -1

    @test Merzbild.vx_sign(4) == 1
    @test Merzbild.vy_sign(4) == 1
    @test Merzbild.vz_sign(4) == -1

    @test Merzbild.vx_sign(5) == -1
    @test Merzbild.vy_sign(5) == -1
    @test Merzbild.vz_sign(5) == 1
    
    @test Merzbild.vx_sign(6) == 1
    @test Merzbild.vy_sign(6) == -1
    @test Merzbild.vz_sign(6) == 1

    @test Merzbild.vx_sign(7) == -1
    @test Merzbild.vy_sign(7) == 1
    @test Merzbild.vz_sign(7) == 1

    @test Merzbild.vx_sign(8) == 1
    @test Merzbild.vy_sign(8) == 1
    @test Merzbild.vz_sign(8) == 1

    @test Merzbild.check_speed_bound(-2.0, 1.0, 3.0) == 1.0
    @test Merzbild.check_speed_bound(5.0, 1.0, 3.0) == 3.0
    @test Merzbild.check_speed_bound(0.0, -1.0, 3.0) == 0.0

    @test Merzbild.check_speed_bound(-20.0, -4.0, 8.0, 0.5) == -2.0
    @test Merzbild.check_speed_bound(20.0, -4.0, 8.0, 0.5) == 4.0
    @test Merzbild.check_speed_bound(-6.0, -4.0, 8.0, 0.5) == -3.0
end