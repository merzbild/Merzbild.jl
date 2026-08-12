@testset "rebalancing" begin
    @testset "ncoll based load balancing init" begin
        lb = LoadBalancerNcoll(100, 4)
        @test lb.n_cells == 100
        @test lb.n_chunks == 4
        @test lb.n_collisions == zeros(Float64, 100)
        @test lb.n_coll_total_per_chunk == zeros(Float64, 4)
        @test lb.n_coll_total == 0.0
        @test lb.chunked_indices == [1:25, 26:50, 51:75, 76:100]

        @test_throws ArgumentError lb2 = LoadBalancerNcoll(3, 4)
    end

    @testset "ncoll based load balancing collision tracking" begin
        lb = LoadBalancerNcoll(6, 4)

        # chunks are 1:3, 4:6, 7:8, 9:10
        update_n_collisions!(lb, 1, 3, 1, 2)  # 3/2
        update_n_collisions!(lb, 1, 0, 2, 2)  # 0
        update_n_collisions!(lb, 1, 2, 3, 2)  # 2/2
        update_n_collisions!(lb, 3, 1, 4, 2)  # 1/2
        update_n_collisions!(lb, 3, 6, 4, 2)  # 6/2

        @test lb.n_coll_total == 0.0 # not set yet
        @test lb.n_collisions == [1.5, 0.0, 1.0, 3.5, 0.0, 0.0]
        @test lb.n_coll_total_per_chunk == [2.5, 0.0, 3.5, 0.0]
        @test sum(lb.n_coll_total_per_chunk) == sum(lb.n_collisions)
        lb.n_coll_total = sum(lb.n_collisions) # set by hand

        reset_lb!(lb)
        @test lb.n_coll_total == 0.0
        @test lb.n_collisions == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    end

    @testset "ncoll based load balancing re-balancing: n_cells=10, n_chunks=4" begin
        lb = LoadBalancerNcoll(10, 4)

        lb.n_collisions = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0]  # 56 collisions

        # 56 collisions in total: 56/4 = 14 per chunk
        # [[1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0], [8.0, 9.0], [11.0]]
        # sums: [15.0, 13.0, 17.0, 11.0]
        rebalance_lb!(lb)
        @test lb.n_coll_total == sum(lb.n_coll_total_per_chunk)
        @test lb.chunked_indices == [1:5, 6:7, 8:9, 10:10]
    end

    @testset "ncoll based load balancing re-balancing: n_cells=10, n_chunks=2" begin
        lb = LoadBalancerNcoll(10, 2)

        lb.n_collisions = [2.0, 3.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # 11 collisions

        # 11 collisions in total: 11/2 = 5.5 per chunk
        # [[2.0, 3.0], [4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        # sums: [5.0, 6.0]
        rebalance_lb!(lb)
        @test lb.chunked_indices == [1:2, 3:10]
    end

    @testset "ncoll based load balancing re-balancing: n_cells=10, n_chunks=4" begin
        lb = LoadBalancerNcoll(10, 4)

        lb.n_collisions = [2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0]  # 11 collisions

        # 11 collisions in total: 11/4 = 2.75 per chunk
        # [[2.0], [3.0], [4.0], [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        # sums: [2.0], [3.0], [4.0], [2.0]
        rebalance_lb!(lb)
        @test lb.chunked_indices == [1:1, 2:2, 3:3, 4:10]
    end

    @testset "ncoll based load balancing re-balancing: n_cells=10, n_chunks=4" begin
        lb = LoadBalancerNcoll(10, 4)

        lb.n_collisions = [0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 4.0, 0.0, 2.0, 0.0]  # 11 collisions

        # 11 collisions in total: 11/4 = 2.75 per chunk
        # sums: [2.0], [3.0], [4.0], [2.0]
        # a chunk is cut off even if some cells at the end have n_coll == 0
        rebalance_lb!(lb)
        @test lb.chunked_indices == [1:4, 5:5, 6:7, 8:10]
    end

    @testset "ncoll based load balancing re-balancing: n_cells=10, n_chunks=4, fitting all chunks in" begin
        lb = LoadBalancerNcoll(10, 4)

        lb.n_collisions = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0] 

        rebalance_lb!(lb)
        @test lb.chunked_indices == [1:7, 8:8, 9:9, 10:10]
    end

    @testset "ncoll based load balancing re-balancing: n_cells=10, n_chunks=4, fitting all chunks in" begin
        lb = LoadBalancerNcoll(10, 4)

        lb.n_collisions = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0] 

        rebalance_lb!(lb)
        @test lb.chunked_indices == [1:7, 8:8, 9:9, 10:10]
    end

    @testset "ncoll based load balancing re-balancing: n_cells=10, n_chunks=4, fitting all chunks in" begin
        lb = LoadBalancerNcoll(10, 4)

        lb.n_collisions = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 20.0]  

        # avg ncoll == 120/10 = 12.0
        # so cutting off the first 6 cells with ncoll_sum == 0.0 is closer

        rebalance_lb!(lb)
        @test lb.chunked_indices == [1:6, 7:7, 8:9, 10:10]
    end
end