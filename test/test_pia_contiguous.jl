@testset "pia_contiguous" begin
    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    # start off with 1-cell 1-group case
    pia = ParticleIndexerArray(1, 1)  # 1 cell 1 species

    particles = [ParticleVector(10)]

    for i in 1:10
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [10.0, 0.0, 1.0])
    end
    particles[1].nbuffer = 0  # set buffer to 0 manually

    # fix particle indexer manually
    pia.indexer[1,1].n_local = 10
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 10
    pia.indexer[1,1].n_group1 = 10

    pia.indexer[1,1].start2 = -1
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0
    
    pia.n_total[1] = 10

    phys_props = PhysProps(pia)
    n_dens = 55.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 10.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    @test pia.contiguous[1] == true

    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)
    pia.contiguous[1] = false

    @test particles[1].nbuffer == 1

    bufferindex = particles[1].buffer[1]

    squash_pia!(particles, pia)
    @test pia.contiguous[1] = true
    
    @test pia.indexer[1,1].n_group1 == 9
    @test pia.indexer[1,1].end1 == 9 

    @test pia.indexer[1,1].n_group2 == 0

    n_dens = 45.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 9.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    # our free particle is not pointed to by any indices in use
    @test (bufferindex in particles[1].index[1:9]) == false


    # 1-cell 2-group case
    for i in 1:10
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [10.0, 0.0, 1.0])
    end
    particles[1].nbuffer = 0  # set buffer to 0 manually

    # fix particle indexer manually
    pia.indexer[1,1].n_local = 10
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 5
    pia.indexer[1,1].n_group1 = 5

    pia.indexer[1,1].start2 = 6
    pia.indexer[1,1].end2 = 10
    pia.indexer[1,1].n_group2 = 5
    
    pia.n_total[1] = 10

    n_dens = 55.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 10.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    @test pia.contiguous[1] == true

    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end_group2!(particles[1], pia, 1, 1)
    pia.contiguous[1] = false

    @test particles[1].nbuffer == 2

    bufferindex1 = particles[1].buffer[1]
    bufferindex2 = particles[1].buffer[2]

    squash_pia!(particles, pia)
    @test pia.contiguous[1] = true
    
    @test pia.indexer[1,1].n_group1 == 4
    @test pia.indexer[1,1].end1 == 4 

    @test pia.indexer[1,1].n_group2 == 4

    n_dens = 40.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 8.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    # our free particles are not pointed to by any indices in use
    @test (bufferindex1 in particles[1].index[1:8]) == false
    @test (bufferindex2 in particles[1].index[1:8]) == false

    # 3-cell case that looks like this

    particles = [ParticleVector(10)]
    # 111 22 33 1 33
    cell_positions = [1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 1.0, 3.0, 3.0]
    for i in 1:10
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [cell_positions[i], 0.0, 1.0])
    end

    particles[1].nbuffer = 0  # set buffer to 0 manually

    pia = ParticleIndexerArray(3, 1)  # 3 cell 1 species
    # fix particle indexer manually
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 3
    pia.indexer[1,1].n_group1 = 3

    pia.indexer[1,1].start2 = 8
    pia.indexer[1,1].end2 = 8
    pia.indexer[1,1].n_group2 = 1

    pia.indexer[2,1].n_local = 2
    pia.indexer[2,1].start1 = 4
    pia.indexer[2,1].end1 = 5
    pia.indexer[2,1].n_group1 = 2

    pia.indexer[2,1].start2 = -1
    pia.indexer[2,1].end2 = -1
    pia.indexer[2,1].n_group2 = 0

    pia.indexer[3,1].n_local = 4
    pia.indexer[3,1].start1 = 6
    pia.indexer[3,1].end1 = 7
    pia.indexer[3,1].n_group1 = 2

    pia.indexer[3,1].start2 = 9
    pia.indexer[3,1].end2 = 10
    pia.indexer[3,1].n_group2 = 2
    
    pia.n_total[1] = 10

    phys_props = PhysProps(pia)
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 4.0
    @test phys_props.np[2,1] == 2.0
    @test phys_props.np[3,1] == 4.0

    @test phys_props.n[1,1] == 14.0
    @test phys_props.n[2,1] == 9.0
    @test phys_props.n[3,1] == 32.0

    # now we delete some stuff so it looks like
    # 1** 22 3* * 33
    # delete 3 times from cell 1
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)

    # delete 1 time from cell 3
    Merzbild.delete_particle_end_group1!(particles[1], pia, 3, 1)
    pia.contiguous[1] = false

    squash_pia!(particles, pia)
    
    # we get
    # 1 22 3 33
    @test pia.indexer[1,1].n_local == 1
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 1
    @test pia.indexer[1,1].n_group1 == 1
    @test pia.indexer[1,1].n_group2 == 0
    
    @test pia.indexer[2,1].n_local == 2
    @test pia.indexer[2,1].start1 == 2
    @test pia.indexer[2,1].end1 == 3
    @test pia.indexer[2,1].n_group1 == 2
    @test pia.indexer[2,1].n_group2 == 0
    
    @test pia.indexer[3,1].n_local == 3
    @test pia.indexer[3,1].start1 == 4
    @test pia.indexer[3,1].end1 == 4
    @test pia.indexer[3,1].n_group1 == 1
    @test pia.indexer[3,1].start2 == 5
    @test pia.indexer[3,1].end2 == 6
    @test pia.indexer[3,1].n_group2 == 2

    for i in 1:particles[1].nbuffer
        @test (particles[1].buffer[i] in particles[1].index[1:6]) == false
    end

    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 1.0
    @test phys_props.np[2,1] == 2.0
    @test phys_props.np[3,1] == 3.0

    @test phys_props.n[1,1] == 1.0
    @test phys_props.n[2,1] == 9.0
    @test phys_props.n[3,1] == 25.0

    # 3-cell case that looks like this
    # 111 2 333 222
    particles = [ParticleVector(10)]
    cell_positions = [1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0]
    for i in 1:10
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [cell_positions[i], 0.0, 1.0])
    end

    particles[1].nbuffer = 0  # set buffer to 0 manually

    pia = ParticleIndexerArray(3, 1)  # 3 cell 1 species
    # fix particle indexer manually
    pia.indexer[1,1].n_local = 3
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 3
    pia.indexer[1,1].n_group1 = 3

    pia.indexer[1,1].start2 = -1
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0

    pia.indexer[2,1].n_local = 4
    pia.indexer[2,1].start1 = 4
    pia.indexer[2,1].end1 = 4
    pia.indexer[2,1].n_group1 = 1

    pia.indexer[2,1].start2 = 8
    pia.indexer[2,1].end2 = 10
    pia.indexer[2,1].n_group2 = 3

    pia.indexer[3,1].n_local = 3
    pia.indexer[3,1].start1 = 5
    pia.indexer[3,1].end1 = 7
    pia.indexer[3,1].n_group1 = 3

    pia.indexer[3,1].start2 = -1
    pia.indexer[3,1].end2 = -1
    pia.indexer[3,1].n_group2 = 0
    
    pia.n_total[1] = 10

    phys_props = PhysProps(pia)
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 3.0
    @test phys_props.np[2,1] == 4.0
    @test phys_props.np[3,1] == 3.0

    @test phys_props.n[1,1] == 6.0
    @test phys_props.n[2,1] == 31.0
    @test phys_props.n[3,1] == 18.0


    # now we delete some stuff so it looks like
    # 11* * 3** 222
    # should still work, although weird that cell 2 has only particles in group2
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 3, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 3, 1)
    Merzbild.delete_particle_end_group1!(particles[1], pia, 2, 1)
    pia.contiguous[1] = false

    squash_pia!(particles, pia)
    # should get
    # 11 3 222
    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].n_group2 == 0
    
    @test pia.indexer[2,1].n_local == 3
    # @test pia.indexer[2,1].start1 == 
    # @test pia.indexer[2,1].end1 == 3
    @test pia.indexer[2,1].n_group1 == 0
    @test pia.indexer[2,1].start2 == 4
    @test pia.indexer[2,1].end2 == 6
    @test pia.indexer[2,1].n_group2 == 3
    
    @test pia.indexer[3,1].n_local == 1
    @test pia.indexer[3,1].start1 == 3
    @test pia.indexer[3,1].end1 == 3
    @test pia.indexer[3,1].n_group1 == 1
    @test pia.indexer[3,1].n_group2 == 0

    # check that there are no conflicts between particles in the buffer and particles
    # actually being used
    flag = false
    for i in 1:6
        if particles[1].index[i] in particles[1].buffer[1:particles[1].nbuffer]
            flag = true
        end
    end
    @test flag == false

    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 2.0
    @test phys_props.np[2,1] == 3.0
    @test phys_props.np[3,1] == 1.0

    @test phys_props.n[1,1] == 3.0
    @test phys_props.n[2,1] == 27.0
    @test phys_props.n[3,1] == 5.0

    # sample 1-D, merge per cell, compute physical props
    particles = [ParticleVector(700)]
    Fnum = 1e10
    pia = ParticleIndexerArray(3, 1)  # 3 cell 1 species
    # sample 100, 200, 300
    xcells = [0.0, 0.5, 1.0, 1.5]
    nps = [100, 200, 300]
    extrabuffer_length = 700 - sum(nps)
    Ts = [1000.0, 750.0, 600.0]

    # reset RNG
    rng = StableRNG(seed)

    for i in 1:3
        sample_particles_equal_weight!(rng, particles[1], pia, i, 1,
                                       nps[i], species_data[1].mass,
                                       Ts[i], Fnum, xcells[i], xcells[i+1], 0.0, 1.0, 0.0, 1.0)
    end
    phys_props = PhysProps(pia)
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.np[1,1] == 100.0
    @test phys_props.np[2,1] == 200.0
    @test phys_props.np[3,1] == 300.0

    @test abs(phys_props.T[1,1] - Ts[1])/Ts[1] < 8.5e-2
    @test abs(phys_props.T[2,1] - Ts[2])/Ts[2] < 1.5e-2
    @test abs(phys_props.T[3,1] - Ts[3])/Ts[3] < 6.5e-2

    T_computed = copy(phys_props.T[:,1])

    octree = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)
    merge_octree_N2_based!(rng, octree, particles[1], pia, 1, 1, 16)
    @test pia.contiguous[1] == false
    merge_octree_N2_based!(rng, octree, particles[1], pia, 2, 1, 16)
    @test pia.contiguous[1] == false
    merge_octree_N2_based!(rng, octree, particles[1], pia, 3, 1, 16)
    @test pia.contiguous[1] == false

    squash_pia!(particles, pia)
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.np[1,1] < 30.0
    @test phys_props.np[2,1] < 30.0
    @test phys_props.np[3,1] < 30.0

    @test abs(phys_props.n[1,1] - Fnum*nps[1])/(Fnum*nps[1]) < 1e-15
    @test abs(phys_props.n[2,1] - Fnum*nps[2])/(Fnum*nps[2]) < 1e-15
    @test abs(phys_props.n[3,1] - Fnum*nps[3])/(Fnum*nps[3]) < 1e-15

    @test abs(phys_props.T[1,1] - T_computed[1])/T_computed[1] < 1e-14
    @test abs(phys_props.T[2,1] - T_computed[2])/T_computed[2] < 1e-14
    @test abs(phys_props.T[3,1] - T_computed[3])/T_computed[3] < 1e-14

    # check that buffer length is equal to the number of unused particles
    # plus whatever buffer we had at the start
    @test particles[1].nbuffer == Int64(sum(nps) - sum(phys_props.np[:,1])) + extrabuffer_length

    # now we do the same, but we split particles of cell 3 into 2 groups
    # sample 1-D, merge per cell, compute physical props
    particles = [ParticleVector(700)]
    Fnum = 1e10
    pia = ParticleIndexerArray(3, 1)  # 3 cell 1 species
    # sample 100, 200, 300
    xcells = [0.0, 0.5, 1.0, 1.5]
    nps = [100, 200, 300]
    extrabuffer_length = 700 - sum(nps)
    Ts = [1000.0, 750.0, 600.0]

    # reset RNG
    rng = StableRNG(seed)

    for i in 1:3
        sample_particles_equal_weight!(rng, particles[1], pia, i, 1,
                                       nps[i], species_data[1].mass,
                                       Ts[i], Fnum, xcells[i], xcells[i+1], 0.0, 1.0, 0.0, 1.0)
    end
    # was 301,...,600
    pia.indexer[3,1].end1 = 304
    pia.indexer[3,1].n_group1 = 4

    pia.indexer[3,1].start2 = 305
    pia.indexer[3,1].end2 = 600
    pia.indexer[3,1].n_group2 = 296
    phys_props = PhysProps(pia)
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.np[1,1] == 100.0
    @test phys_props.np[2,1] == 200.0
    @test phys_props.np[3,1] == 300.0

    @test abs(phys_props.T[1,1] - Ts[1])/Ts[1] < 8.5e-2
    @test abs(phys_props.T[2,1] - Ts[2])/Ts[2] < 1.5e-2
    @test abs(phys_props.T[3,1] - Ts[3])/Ts[3] < 6.5e-2

    merge_octree_N2_based!(rng, octree, particles[1], pia, 3, 1, 16)
    @test pia.contiguous[1] == true
    merge_octree_N2_based!(rng, octree, particles[1], pia, 2, 1, 16)
    @test pia.contiguous[1] == false
    merge_octree_N2_based!(rng, octree, particles[1], pia, 1, 1, 16)
    @test pia.contiguous[1] == false

    squash_pia!(particles, pia)
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.np[1,1] < 30.0
    @test phys_props.np[2,1] < 30.0
    @test phys_props.np[3,1] < 30.0

    @test abs(phys_props.n[1,1] - Fnum*nps[1])/(Fnum*nps[1]) < 1e-15
    @test abs(phys_props.n[2,1] - Fnum*nps[2])/(Fnum*nps[2]) < 1e-15
    @test abs(phys_props.n[3,1] - Fnum*nps[3])/(Fnum*nps[3]) < 1e-15

    @test abs(phys_props.T[1,1] - T_computed[1])/T_computed[1] < 1e-14
    @test abs(phys_props.T[2,1] - T_computed[2])/T_computed[2] < 1e-14
    @test abs(phys_props.T[3,1] - T_computed[3])/T_computed[3] < 1e-14


    # 3-cell case that looks like this
    # 111 2 333 222
    particles = [ParticleVector(10)]
    cell_positions = [1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0]
    for i in 1:10
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [cell_positions[i], 0.0, 1.0])
    end

    particles[1].nbuffer = 0  # set buffer to 0 manually

    pia = ParticleIndexerArray(3, 1)  # 3 cell 1 species
    # fix particle indexer manually
    pia.indexer[1,1].n_local = 3
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 3
    pia.indexer[1,1].n_group1 = 3

    pia.indexer[1,1].start2 = -1
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0

    pia.indexer[2,1].n_local = 4
    pia.indexer[2,1].start1 = 4
    pia.indexer[2,1].end1 = 4
    pia.indexer[2,1].n_group1 = 1

    pia.indexer[2,1].start2 = 8
    pia.indexer[2,1].end2 = 10
    pia.indexer[2,1].n_group2 = 3

    pia.indexer[3,1].n_local = 3
    pia.indexer[3,1].start1 = 5
    pia.indexer[3,1].end1 = 7
    pia.indexer[3,1].n_group1 = 3

    pia.indexer[3,1].start2 = -1
    pia.indexer[3,1].end2 = -1
    pia.indexer[3,1].n_group2 = 0
    
    pia.n_total[1] = 10

    phys_props = PhysProps(pia)
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 3.0
    @test phys_props.np[2,1] == 4.0
    @test phys_props.np[3,1] == 3.0

    @test phys_props.n[1,1] == 6.0
    @test phys_props.n[2,1] == 31.0
    @test phys_props.n[3,1] == 18.0


    # now we delete some stuff so it looks like
    # *** 2 333 222
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    pia.contiguous[1] = false

    squash_pia!(particles, pia)
    # should get
    # 2 333 222
    @test pia.indexer[1,1].n_local == 0
    @test pia.indexer[1,1].start1 == 0
    @test pia.indexer[1,1].end1 == -1
    @test pia.indexer[1,1].n_group1 == 0
    @test pia.indexer[1,1].n_group2 == 0
    
    @test pia.indexer[2,1].n_local == 4
    @test pia.indexer[2,1].n_group1 == 1
    @test pia.indexer[2,1].start1 == 1
    @test pia.indexer[2,1].end1 == 1
    @test pia.indexer[2,1].start2 == 5
    @test pia.indexer[2,1].end2 == 7
    @test pia.indexer[2,1].n_group2 == 3
    
    @test pia.indexer[3,1].n_local == 3
    @test pia.indexer[3,1].start1 == 2
    @test pia.indexer[3,1].end1 == 4
    @test pia.indexer[3,1].n_group1 == 3
    @test pia.indexer[3,1].n_group2 == 0

    # more tests
    # 3-cell case that looks like this
    # 111 2 333 222
    particles = [ParticleVector(10)]
    cell_positions = [1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0]
    for i in 1:10
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [cell_positions[i], 0.0, 1.0])
    end

    particles[1].nbuffer = 0  # set buffer to 0 manually

    pia = ParticleIndexerArray(3, 1)  # 3 cell 1 species
    # fix particle indexer manually
    pia.indexer[1,1].n_local = 3
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 3
    pia.indexer[1,1].n_group1 = 3

    pia.indexer[1,1].start2 = -1
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0

    pia.indexer[2,1].n_local = 4
    pia.indexer[2,1].start1 = 4
    pia.indexer[2,1].end1 = 4
    pia.indexer[2,1].n_group1 = 1

    pia.indexer[2,1].start2 = 8
    pia.indexer[2,1].end2 = 10
    pia.indexer[2,1].n_group2 = 3

    pia.indexer[3,1].n_local = 3
    pia.indexer[3,1].start1 = 5
    pia.indexer[3,1].end1 = 7
    pia.indexer[3,1].n_group1 = 3

    pia.indexer[3,1].start2 = -1
    pia.indexer[3,1].end2 = -1
    pia.indexer[3,1].n_group2 = 0
    
    pia.n_total[1] = 10

    phys_props = PhysProps(pia)
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 3.0
    @test phys_props.np[2,1] == 4.0
    @test phys_props.np[3,1] == 3.0

    @test phys_props.n[1,1] == 6.0
    @test phys_props.n[2,1] == 31.0
    @test phys_props.n[3,1] == 18.0


    # now we delete some stuff so it looks like
    # 11* 2 *** 222
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 3, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 3, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 3, 1)
    pia.contiguous[1] = false

    squash_pia!(particles, pia)
    # should get
    # 11 2 222
    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].n_group2 == 0
    
    @test pia.indexer[2,1].n_local == 4
    @test pia.indexer[2,1].n_group1 == 1
    @test pia.indexer[2,1].start1 == 3
    @test pia.indexer[2,1].end1 == 3
    @test pia.indexer[2,1].start2 == 4
    @test pia.indexer[2,1].end2 == 6
    @test pia.indexer[2,1].n_group2 == 3
    
    @test pia.indexer[3,1].n_local == 0
    @test pia.indexer[3,1].start1 == 0
    @test pia.indexer[3,1].end1 == -1
    @test pia.indexer[3,1].n_group1 == 0
    @test pia.indexer[3,1].n_group2 == 0
end