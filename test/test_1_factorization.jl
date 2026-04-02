@testset "test graph factorization" begin

    fac2 = generate_1_factorization(2)

    @test length(fac2) == 1
    @test length(fac2[1]) == 1
    @test ((fac2[1][1] == (1,2)) || (fac2[1][1] == (2,1))) == true

    for N_elems in [4,8,16,32]
        fac_long = generate_1_factorization(N_elems)

        @test length(fac_long) == N_elems - 1

        # at least at the start we have maximum possible length of vectors
        @test length(fac_long[1]) == N_elems ÷ 2
        @test length(fac_long[2]) == N_elems ÷ 2

        # each vector has unique elements in the tuples
        for row in fac_long
            element_count = zeros(N_elems)
            for tp in row
                element_count[tp[1]] += 1
                element_count[tp[2]] += 1
            end

            @test maximum(element_count) == 1
        end

        # all pairs are there
        for i in 1:N_elems
            for j in 1:N_elems
                if i != j
                    found = false
                    for row in fac_long
                        for tp in row
                            if ((i,j) == tp || (j,i) == tp)
                                found = true
                                break
                            end
                        end
                    end
                    @test found == true
                else 
                    # should not show up!
                    found = false
                    for row in fac_long
                        for tp in row
                            if ((i,i) == tp || (i,i) == tp)
                                found = true
                                break
                            end
                        end
                    end
                    @test found == false
                end
            end
        end
    end

    # and for other counts we also get what we expect
    for N_elems in [5,20,31]
        fac_long = generate_1_factorization(N_elems)

        # at least at the start we have maximum possible length of vectors
        @test length(fac_long[1]) == N_elems ÷ 2
        @test length(fac_long[2]) == N_elems ÷ 2

        # each vector has unique elements in the tuples
        for row in fac_long
            element_count = zeros(N_elems)
            for tp in row
                element_count[tp[1]] += 1
                element_count[tp[2]] += 1
            end

            @test maximum(element_count) == 1
        end

        # all pairs are there
        for i in 1:N_elems
            for j in 1:N_elems
                if i != j
                    found = false
                    for row in fac_long
                        for tp in row
                            if ((i,j) == tp || (j,i) == tp)
                                found = true
                                break
                            end
                        end
                    end
                    @test found == true
                else 
                    # should not show up!
                    found = false
                    for row in fac_long
                        for tp in row
                            if ((i,i) == tp || (i,i) == tp)
                                found = true
                                break
                            end
                        end
                    end
                    @test found == false
                end
            end
        end
    end
end