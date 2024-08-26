@testset "misc" begin
    x = [1.0, 2.0, 3.5]
    val = 1.1
    @test Merzbild.binary_search(x, val) == 1

    val = 2.0
    @test Merzbild.binary_search(x, val) == 2

    val = 5.0
    @test Merzbild.binary_search(x, val) == 0

    val = 0.9
    @test Merzbild.binary_search(x, val) == -1

    x = [1, 3, 4, 5]
    weights = ones(length(x))
    mymedian = Merzbild.weighted_percentile_interpolated(x, weights, 0.5)
    @test abs(mymedian - 3.0) < eps()

    x = [0.5, 1, 2, 4, 5]
    weights = ones(length(x))
    mymedian = Merzbild.weighted_percentile_interpolated(x, weights, 0.5)
    @test abs(mymedian - 1.5) < eps()
end