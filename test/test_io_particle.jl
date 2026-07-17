@testset "I/O of particles" begin
    # various basic I/O routines
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, ["Ar", "He"])
    
    sol_path = joinpath(@__DIR__, "data", "tmp_1species.nc")
    
    pia = ParticleIndexerArray(3,2)

    particles = [ParticleVector(5), ParticleVector(6)]

    #    cell 1,   2,   2    3     2
    x_pos1 = [0.0, 1.0, 3.0, 10.0, 4.0]
    for i in 1:5
        particles[1][i] = Particle(i*1.0, [i + 2.0, -i, i*2.0], [x_pos1[i], x_pos1[i] * -1, x_pos1[i] + 3.0])
    end

    #    cell 1,   2,   2    3     2    3
    x_pos2 = [1.0, 2.0, 4.0, 11.0, 3.0, 12.0]
    for i in 1:6
        particles[2][i] = Particle(3*i+10.0, [-i - 4.0, i, i*4.0 + 5.0], [x_pos2[i], x_pos2[i] * 2 + 1.0, -x_pos2[i] - 3.0])
    end

    pia.n_total[1] = 5
    pia.n_total[2] = 6

    pia.indexer[1,1].n_local = 1
    pia.indexer[1,1].n_group1 = 1
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 1

    pia.indexer[2,1].n_local = 3
    pia.indexer[2,1].n_group1 = 2
    pia.indexer[2,1].start1 = 2
    pia.indexer[2,1].end1 = 3

    pia.indexer[2,1].n_group2 = 1
    pia.indexer[2,1].start2 = 5
    pia.indexer[2,1].end2 = 5

    pia.indexer[3,1].n_local = 1
    pia.indexer[3,1].n_group1 = 1
    pia.indexer[3,1].start1 = 4
    pia.indexer[3,1].end1 = 4

    # particles of species 2
    pia.indexer[1,2].n_local = 1
    pia.indexer[1,2].n_group1 = 1
    pia.indexer[1,2].start1 = 1
    pia.indexer[1,2].end1 = 1

    pia.indexer[2,2].n_local = 3
    pia.indexer[2,2].n_group1 = 2
    pia.indexer[2,2].start1 = 2
    pia.indexer[2,2].end1 = 3

    pia.indexer[2,2].n_group2 = 1
    pia.indexer[2,2].start2 = 5
    pia.indexer[2,2].end2 = 5

    pia.indexer[3,2].n_local = 2
    pia.indexer[3,2].n_group1 = 1
    pia.indexer[3,2].start1 = 4
    pia.indexer[3,2].end1 = 4
    pia.indexer[3,2].n_group2 = 1
    pia.indexer[3,2].start2 = 6
    pia.indexer[3,2].end2 = 6

    write_netcdf(sol_path, particles, pia, species_data, [1])

    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in ["w_Ar", "v_Ar", "x_Ar", "cell_Ar"]
        @test name in varnames
    end

    @test size(props_read["x_Ar"]) == (5, 3)
    @test size(props_read["v_Ar"]) == (5, 3)
    @test size(props_read["w_Ar"]) == (5,)
    @test size(props_read["cell_Ar"]) == (5,)
    @test props_read.attrib["species_names"] == "Ar"
    @test props_read.attrib["Particle x-position dimension"] == "3"

    # we write cell-by-cell, so some re-ordering happens
    @test props_read["w_Ar"][:] == [1.0, 2.0, 3.0, 5.0, 4.0]

    @test props_read["x_Ar"][:, 1] == [0.0, 1.0, 3.0, 4.0, 10.0]
    @test props_read["x_Ar"][:, 2] == -1.0 .* [0.0, 1.0, 3.0, 4.0, 10.0]
    @test props_read["x_Ar"][:, 3] == [0.0, 1.0, 3.0, 4.0, 10.0] .+ 3.0
 
    @test props_read["v_Ar"][:, 1] == [1.0, 2.0, 3.0, 5.0, 4.0] .+ 2.0
    @test props_read["v_Ar"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0] * -1.0
    @test props_read["v_Ar"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0] * 2.0

    @test props_read["cell_Ar"][:, 1] == [1, 2, 2, 2, 3]

    close(props_read)
    rm(sol_path)

    sol_path = joinpath(@__DIR__, "data", "tmp_all_species.nc")
    write_netcdf(sol_path, particles, pia, species_data, [1,2])
    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in ["w_Ar", "v_Ar", "x_Ar", "cell_Ar", "w_He", "v_He", "x_He", "cell_He"]
        @test name in varnames
    end

    @test props_read.attrib["species_names"] == "Ar,He"
    @test props_read.attrib["Particle x-position dimension"] == "3"

    @test props_read["w_Ar"][:] == [1.0, 2.0, 3.0, 5.0, 4.0]

    @test props_read["x_Ar"][:, 1] == [0.0, 1.0, 3.0, 4.0, 10.0]
    @test props_read["x_Ar"][:, 2] == -1.0 .* [0.0, 1.0, 3.0, 4.0, 10.0]
    @test props_read["x_Ar"][:, 3] == [0.0, 1.0, 3.0, 4.0, 10.0] .+ 3.0
 
    @test props_read["v_Ar"][:, 1] == [1.0, 2.0, 3.0, 5.0, 4.0] .+ 2.0
    @test props_read["v_Ar"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0] * -1.0
    @test props_read["v_Ar"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0] * 2.0

    @test props_read["cell_Ar"][:, 1] == [1, 2, 2, 2, 3]

    @test props_read["w_He"][:] == 3.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .+ 10.0

    @test props_read["x_He"][:, 1] == [1.0, 2.0, 4.0, 3.0, 11.0, 12.0]
    @test props_read["x_He"][:, 2] == 1.0 .+ 2.0 .* [1.0, 2.0, 4.0, 3.0, 11.0, 12.0]
    @test props_read["x_He"][:, 3] == -1.0 * [1.0, 2.0, 4.0, 3.0, 11.0, 12.0] .- 3.0
 
    @test props_read["v_He"][:, 1] == -1.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .- 4.0
    @test props_read["v_He"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0]
    @test props_read["v_He"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] * 4.0 .+ 5.0

    @test props_read["cell_He"][:, 1] == [1, 2, 2, 2, 3, 3]

    close(props_read)
    rm(sol_path)

    # test 2D particles
    particles = [ParticleVector{2}(5), ParticleVector{2}(6)]

    #    cell 1,   2,   2    3     2
    x_pos1 = [0.0, 1.0, 3.0, 10.0, 4.0]
    for i in 1:5
        particles[1][i] = Particle{2}(i*1.0, SVector{3,Float64}(i + 2.0, -i, i*2.0), SVector{2,Float64}(x_pos1[i], x_pos1[i] * -1))
    end

    #    cell 1,   2,   2    3     2    3
    x_pos2 = [1.0, 2.0, 4.0, 11.0, 3.0, 12.0]
    for i in 1:6
        particles[2][i] = Particle{2}(3*i+10.0, SVector{3,Float64}(-i - 4.0, i, i*4.0 + 5.0), SVector{2,Float64}(x_pos2[i], x_pos2[i] * 2 + 1.0))
    end


    sol_path = joinpath(@__DIR__, "data", "tmp_all_species.nc")
    write_netcdf(sol_path, particles, pia, species_data, [1,2])
    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in ["w_Ar", "v_Ar", "x_Ar", "cell_Ar", "w_He", "v_He", "x_He", "cell_He"]
        @test name in varnames
    end

    @test props_read.attrib["species_names"] == "Ar,He"
    @test props_read.attrib["Particle x-position dimension"] == "2"

    @test props_read["w_Ar"][:] == [1.0, 2.0, 3.0, 5.0, 4.0]

    @test props_read["x_Ar"][:, 1] == [0.0, 1.0, 3.0, 4.0, 10.0]
    @test props_read["x_Ar"][:, 2] == -1.0 .* [0.0, 1.0, 3.0, 4.0, 10.0]
 
    @test props_read["v_Ar"][:, 1] == [1.0, 2.0, 3.0, 5.0, 4.0] .+ 2.0
    @test props_read["v_Ar"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0] * -1.0
    @test props_read["v_Ar"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0] * 2.0

    @test props_read["cell_Ar"][:, 1] == [1, 2, 2, 2, 3]

    @test props_read["w_He"][:] == 3.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .+ 10.0

    @test props_read["x_He"][:, 1] == [1.0, 2.0, 4.0, 3.0, 11.0, 12.0]
    @test props_read["x_He"][:, 2] == 1.0 .+ 2.0 .* [1.0, 2.0, 4.0, 3.0, 11.0, 12.0]
 
    @test props_read["v_He"][:, 1] == -1.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .- 4.0
    @test props_read["v_He"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0]
    @test props_read["v_He"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] * 4.0 .+ 5.0

    @test props_read["cell_He"][:, 1] == [1, 2, 2, 2, 3, 3]

    close(props_read)
    rm(sol_path)

    # test 1D particles
    particles = [ParticleVector{1}(5), ParticleVector{1}(6)]

    #    cell 1,   2,   2    3     2
    x_pos1 = [0.0, 1.0, 3.0, 10.0, 4.0]
    for i in 1:5
        particles[1][i] = Particle{1}(i*1.0, SVector{3,Float64}(i + 2.0, -i, i*2.0), SVector{1,Float64}(x_pos1[i]))
    end

    #    cell 1,   2,   2    3     2    3
    x_pos2 = [1.0, 2.0, 4.0, 11.0, 3.0, 12.0]
    for i in 1:6
        particles[2][i] = Particle{1}(3*i+10.0, SVector{3,Float64}(-i - 4.0, i, i*4.0 + 5.0), SVector{1,Float64}(x_pos2[i]))
    end

    sol_path = joinpath(@__DIR__, "data", "tmp_all_species.nc")
    write_netcdf(sol_path, particles, pia, species_data, [1,2])
    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in ["w_Ar", "v_Ar", "x_Ar", "cell_Ar", "w_He", "v_He", "x_He", "cell_He"]
        @test name in varnames
    end

    @test props_read.attrib["species_names"] == "Ar,He"
    @test props_read.attrib["Particle x-position dimension"] == "1"

    @test props_read["w_Ar"][:] == [1.0, 2.0, 3.0, 5.0, 4.0]

    @test props_read["x_Ar"][:, 1] == [0.0, 1.0, 3.0, 4.0, 10.0]
 
    @test props_read["v_Ar"][:, 1] == [1.0, 2.0, 3.0, 5.0, 4.0] .+ 2.0
    @test props_read["v_Ar"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0] * -1.0
    @test props_read["v_Ar"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0] * 2.0

    @test props_read["cell_Ar"][:, 1] == [1, 2, 2, 2, 3]

    @test props_read["w_He"][:] == 3.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .+ 10.0

    @test props_read["x_He"][:, 1] == [1.0, 2.0, 4.0, 3.0, 11.0, 12.0]
 
    @test props_read["v_He"][:, 1] == -1.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .- 4.0
    @test props_read["v_He"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0]
    @test props_read["v_He"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] * 4.0 .+ 5.0

    @test props_read["cell_He"][:, 1] == [1, 2, 2, 2, 3, 3]

    close(props_read)
    rm(sol_path)

    # test 0D particles
    particles = [ParticleVector{0}(5), ParticleVector{0}(6)]

    #    cell 1,   2,   2    3     2
    x_pos1 = [0.0, 1.0, 3.0, 10.0, 4.0]
    for i in 1:5
        particles[1][i] = Particle{0}(i*1.0, SVector{3,Float64}(i + 2.0, -i, i*2.0), SVector{0,Float64}())
    end

    #    cell 1,   2,   2    3     2    3
    x_pos2 = [1.0, 2.0, 4.0, 11.0, 3.0, 12.0]
    for i in 1:6
        particles[2][i] = Particle{0}(3*i+10.0, SVector{3,Float64}(-i - 4.0, i, i*4.0 + 5.0), SVector{0,Float64}())
    end

    sol_path = joinpath(@__DIR__, "data", "tmp_all_species.nc")
    write_netcdf(sol_path, particles, pia, species_data, [1,2])
    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in ["w_Ar", "v_Ar", "w_He", "v_He", ]
        @test name in varnames
    end

    for name in ["x_Ar", "cell_Ar", "x_He", "cell_He"]
        @test !(name in varnames)
    end

    @test props_read.attrib["species_names"] == "Ar,He"
    @test props_read.attrib["Particle x-position dimension"] == "0"

    @test props_read["w_Ar"][:] == [1.0, 2.0, 3.0, 5.0, 4.0]
 
    @test props_read["v_Ar"][:, 1] == [1.0, 2.0, 3.0, 5.0, 4.0] .+ 2.0
    @test props_read["v_Ar"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0] * -1.0
    @test props_read["v_Ar"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0] * 2.0

    @test props_read["w_He"][:] == 3.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .+ 10.0
 
    @test props_read["v_He"][:, 1] == -1.0 * [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] .- 4.0
    @test props_read["v_He"][:, 2] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0]
    @test props_read["v_He"][:, 3] == [1.0, 2.0, 3.0, 5.0, 4.0, 6.0] * 4.0 .+ 5.0


    close(props_read)
    rm(sol_path)
end