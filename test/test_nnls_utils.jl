@testset "nnls_merging utils" begin
    base_mim = Merzbild.base_multi_index_moments()
    ref_mom_indices = [[0, 0, 0],
                       [1, 0, 0], [0, 1, 0], [0, 0, 1],
                       [2, 0, 0], [0, 2, 0], [0, 0, 2]]

    @test length(base_mim) == length(ref_mom_indices)
    for ref_index in ref_mom_indices
        @test ref_index in base_mim
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

    mim = Merzbild.compute_multi_index_moments(3)

    for x in mim
        @test sum(x) == 3
    end
    @test [3,0,0] in mim
    @test [0,3,0] in mim
    @test [0,0,3] in mim
    @test [2,1,0] in mim
    @test [0,1,2] in mim
    @test [1,0,2] in mim
    @test [2,0,1] in mim
    @test [0,2,1] in mim
    @test [1,2,0] in mim
    @test [1,1,1] in mim
    @test length(mim) == 10

    matrix = [3.0 5.0 7.0; 4.0 12.0 24.0]
    column_norms = [0.0, 0.0, 0.0]

    Merzbild.scale_columns!(matrix, column_norms)

    for i in 1:3
        @test abs(norm(matrix[:,i], 2) - 1.0) < 2*eps()
    end
    @test abs(column_norms[1] - 1.0/5.0) < 2*eps()
    @test abs(column_norms[2] - 1.0/13.0) < 2*eps()
    @test abs(column_norms[3] - 1.0/25.0) < 2*eps()
end