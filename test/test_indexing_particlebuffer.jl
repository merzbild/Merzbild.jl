@testset "particle buffer" begin
    # the tests of indexing ParticleVectors and the buffer part of the ParticleVectors
    # are not that easy to split nicely, so some duplication of tests is possible

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")
    pia = ParticleIndexerArray(0)
    
    particles = [ParticleVector(10)]

    seed = 1234
    rng = StableRNG(seed)

    

    # test that for non-initialized particles the buffer is correct
    @test particles[1].buffer == [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    @test particles[1].nbuffer == 10

    nparticles = 6

    # mass = 1.0, fnum = 1e10
    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1,
                                   nparticles, species_data[1].mass,
                                   237.0, 1e10, 0, 1.0, 0, 1.0, 0.0, 1.0;
                                   distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)

    for i in 1:6
        @test particles[1][i].w == 1e10
    end
    for i in 7:10
        @test isdefined(particles[1], i) == false
    end

    @test particles[1].nbuffer == 10 - nparticles
    @test particles[1].buffer == [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

    # add 4 new particles at the end
    resize!(particles[1], 14)
    @test particles[1].index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    @test particles[1].nbuffer == 10 - nparticles + 4

    # new indices are added in reverse at the start of the buffer array (older particles get used up first)
    @test particles[1].buffer == [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]


    # mass = 1.0, fnum = 16
    nparticles2 = 3
    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1,
                                   nparticles2, 1.0, 237.0, 1e6, 0, 1.0, 0, 1.0, 0.0, 1.0;
                                   distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)
    
    @test particles[1].nbuffer == 10 - nparticles + 4 - nparticles2
    lenbuf = particles[1].nbuffer

    for i in 1:6
        @test particles[1][i].w == 1e10
    end
    for i in 7:9
        @test particles[1][i].w == 1e6
    end
    for i in 10:14
        isdefined(particles[1], i) == false
    end

    # fix particle indexer manually
    pia.indexer[1,1].n_local = 9
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 6
    pia.indexer[1,1].n_group1 = 6

    pia.indexer[1,1].start2 = 7
    pia.indexer[1,1].end2 = 9
    pia.indexer[1,1].n_group2 = 3
    
    pia.n_total[1] = 9

    n_dens = 6e10 + 3e6

    phys_props = PhysProps(1, 1, [], Tref=1)
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 9.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    # delete particle from end of group1
    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)

    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 5
    @test particles[1].nbuffer == lenbuf + 1
    # our buffer looked like this: [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    # and we had 5 elements there, [14, 13, 12, 11, 10]
    # now we removed particle with index 6
    @test particles[1].buffer == [14, 13, 12, 11, 10, 6, 8, 7, 6, 5, 4, 3, 2, 1]
    @test particles[1].index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    n_dens = 5e10 + 3e6
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 8.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    # delete particle from end of group2
    Merzbild.delete_particle_end_group2!(particles[1], pia, 1, 1)

    @test particles[1].nbuffer == lenbuf + 2
    @test particles[1].buffer == [14, 13, 12, 11, 10, 6, 9, 7, 6, 5, 4, 3, 2, 1]

    n_dens = 5e10 + 2e6
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 7.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    Merzbild.delete_particle_end_group2!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end_group2!(particles[1], pia, 1, 1)

    @test particles[1].buffer == [14, 13, 12, 11, 10, 6, 9, 8, 7, 5, 4, 3, 2, 1]
    @test particles[1].nbuffer == lenbuf + 4

    n_dens = 5e10
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 5.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15
    @test pia.indexer[1,1].start2 == -1
    @test pia.indexer[1,1].n_group2 == 0

    # now we try deleting again and nothing changes
    Merzbild.delete_particle_end_group2!(particles[1], pia, 1, 1)

    @test particles[1].buffer == [14, 13, 12, 11, 10, 6, 9, 8, 7, 5, 4, 3, 2, 1]
    @test particles[1].nbuffer == lenbuf + 4

    n_dens = 5e10
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 5.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15
    @test pia.indexer[1,1].start2 == -1
    @test pia.indexer[1,1].n_group2 == 0


    # now we try deleting from group1
    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)

    @test particles[1].buffer == [14, 13, 12, 11, 10, 6, 9, 8, 7, 5, 4, 3, 2, 1]
    @test particles[1].nbuffer == lenbuf + 8

    n_dens = 1e10
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 1.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 1
    @test pia.indexer[1,1].start2 == -1
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].n_group1 == 1


    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)

    @test particles[1].buffer == [14, 13, 12, 11, 10, 6, 9, 8, 7, 5, 4, 3, 2, 1]
    @test particles[1].nbuffer == 14

    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 0.0
    @test abs(phys_props.n[1,1] - 0.0) < 1e-15
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 0
    @test pia.indexer[1,1].start2 == -1
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].n_group1 == 0


    # nothing changes, since there are no particles to delete
    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)

    @test particles[1].buffer == [14, 13, 12, 11, 10, 6, 9, 8, 7, 5, 4, 3, 2, 1]
    @test particles[1].nbuffer == 14

    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 0.0
    @test abs(phys_props.n[1,1] - 0.0) < 1e-15
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 0
    @test pia.indexer[1,1].start2 == -1
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].n_group1 == 0

    # test other deletion functions, but also with a non-sorted index array
    # reset pia manually
    pia.indexer[1,1].n_local = 7
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 3
    pia.indexer[1,1].n_group1 = 3

    pia.indexer[1,1].start2 = 7
    pia.indexer[1,1].end2 = 10
    pia.indexer[1,1].n_group2 = 4
    
    pia.n_total[1] = 7

    # write the indices in some unsorted order
    particles[1].index = [10, 7, 9, 12, 1, 5, 11, 2, 13, 3, 4, 6, 14, 8]

    # the first 7 indices are the particles not used at the moment
    particles[1].buffer = [12, 1, 5, 4, 6, 14, 8, 10, 7, 9, 11, 2, 13, 3]
    particles[1].nbuffer = 7

    for ind in particles[1].index
        particles[1].particles[ind] = Particle(ind * 1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    end

    n_dens = 55.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 7.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    Merzbild.delete_particle_end_group1!(particles[1], pia, 1, 1)
    @test particles[1].index == [10, 7, 9, 12, 1, 5, 11, 2, 13, 3, 4, 6, 14, 8]
    @test particles[1].buffer == [12, 1, 5, 4, 6, 14, 8, 9, 7, 9, 11, 2, 13, 3]
    @test particles[1].nbuffer == 8

    # deleted particle with index 9, so ndens is decreased by 9.0
    n_dens = 46.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 6.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    Merzbild.delete_particle_end_group2!(particles[1], pia, 1, 1)
    @test particles[1].index == [10, 7, 9, 12, 1, 5, 11, 2, 13, 3, 4, 6, 14, 8]
    @test particles[1].buffer == [12, 1, 5, 4, 6, 14, 8, 9, 3, 9, 11, 2, 13, 3]
    @test particles[1].nbuffer == 9

    # deleted particle with index 3, so ndens is decreased by 3
    n_dens = 43.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 5.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15

    # this deletes particle at the end of group 2
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    
    @test particles[1].index == [10, 7, 9, 12, 1, 5, 11, 2, 13, 3, 4, 6, 14, 8]
    @test particles[1].buffer == [12, 1, 5, 4, 6, 14, 8, 9, 3, 13, 11, 2, 13, 3]
    @test particles[1].nbuffer == 10

    # deleted particle with index 13, so ndens is decreased by 13
    n_dens = 30.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 4.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15


    # our particle indices (the ones in pia) are [1,2], [7,8]
    # let's delete particle 1
    Merzbild.delete_particle!(particles[1], pia, 1, 1, 1)
    
    # swapped particle-to-be-deleted with last particle in group
    @test particles[1].index == [7, 10, 9, 12, 1, 5, 11, 2, 13, 3, 4, 6, 14, 8]
    @test particles[1].buffer == [12, 1, 5, 4, 6, 14, 8, 9, 3, 13, 10, 2, 13, 3]
    @test particles[1].nbuffer == 11

    # deleted particle with index 10, so ndens is decreased by 10
    n_dens = 20.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 3.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15


    # our particle indices (the ones in pia) are [1], [7,8]
    # let's delete particle 7
    Merzbild.delete_particle!(particles[1], pia, 1, 1, 7)
    
    # swapped particle-to-be-deleted with last particle in group
    @test particles[1].index == [7, 10, 9, 12, 1, 5, 2, 11, 13, 3, 4, 6, 14, 8]
    @test particles[1].buffer == [12, 1, 5, 4, 6, 14, 8, 9, 3, 13, 10, 11, 13, 3]
    @test particles[1].nbuffer == 12

    # deleted particle with index 11, so ndens is decreased by 11
    n_dens = 9.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 2.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15


    # our particle indices (the ones in pia) are [1], [7]
    # let's delete particle 7
    Merzbild.delete_particle!(particles[1], pia, 1, 1, 7)
    
    # swapped particle-to-be-deleted with last particle in group
    @test particles[1].index == [7, 10, 9, 12, 1, 5, 2, 11, 13, 3, 4, 6, 14, 8]
    @test particles[1].buffer == [12, 1, 5, 4, 6, 14, 8, 9, 3, 13, 10, 11, 2, 3]
    @test particles[1].nbuffer == 13

    # deleted particle with index 2, so ndens is decreased by 2
    n_dens = 7.0
    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 1.0
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < 1e-15


    # this deletes particle at the end of group 1
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    @test particles[1].index == [7, 10, 9, 12, 1, 5, 2, 11, 13, 3, 4, 6, 14, 8]
    @test particles[1].buffer == [12, 1, 5, 4, 6, 14, 8, 9, 3, 13, 10, 11, 2, 7]
    @test particles[1].nbuffer == 14

    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 0.0
    @test abs(phys_props.n[1,1]) < 1e-15
end