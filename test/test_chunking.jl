@testset "chunk operations and 1-D chunk sampling" begin
    # test that computing properties on chunks works correctly
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(10.0, 5)

    particles = [ParticleVector(4)]
    rng = StableRNG(1234)

    cell_positions = [3.0, 3.0, 3.0, 7.0]
    for i in 1:4
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [cell_positions[i], 0.0, 1.0])
    end

    particles[1].nbuffer = 0  # set buffer to 0 manually

    pia = ParticleIndexerArray(grid.n_cells, 1)
    phys_props = PhysProps(pia)

    # fix particle indexer manually
    pia.indexer[2,1].n_local = 3
    pia.indexer[2,1].start1 = 1
    pia.indexer[2,1].end1 = 3
    pia.indexer[2,1].n_group1 = 3
    pia.indexer[2,1].start2 = 0
    pia.indexer[2,1].end2 = -1
    pia.indexer[2,1].n_group2 = 0

    pia.indexer[3,1].n_local = 0
    pia.indexer[3,1].start1 = 0
    pia.indexer[3,1].end1 = -1
    pia.indexer[3,1].n_group1 = 0
    pia.indexer[3,1].start2 = 0
    pia.indexer[3,1].end2 = -1
    pia.indexer[3,1].n_group2 = 0

    pia.indexer[4,1].n_local = 1
    pia.indexer[4,1].start1 = 4
    pia.indexer[4,1].end1 = 4
    pia.indexer[4,1].n_group1 = 1
    pia.indexer[4,1].start2 = 0
    pia.indexer[4,1].end2 = -1
    pia.indexer[4,1].n_group2 = 0

    for i in [1,3,5]
        pia.indexer[i,1].n_local = 0
        pia.indexer[i,1].start1 = 0
        pia.indexer[i,1].end1 = -1
        pia.indexer[i,1].n_group1 = 0
        pia.indexer[i,1].start2 = 0
        pia.indexer[i,1].end2 = -1
        pia.indexer[i,1].n_group2 = 0
    end

    pia.n_total[1] = 4

    cell_chunk = [1]
    compute_props_sorted!(particles, pia, species_data, phys_props, cell_chunk)
    for i in 1:5
        @test phys_props.np[i,1] == 0.0
    end

    cell_chunk = [1,2,3]
    compute_props_sorted!(particles, pia, species_data, phys_props, cell_chunk)
    for i in [1,3,4,5]
        @test phys_props.np[i,1] == 0.0
    end
    @test phys_props.np[2,1] == 3.0
    @test phys_props.n[2,1] == 6.0

    cell_chunk = [3,4]
    compute_props_sorted!(particles, pia, species_data, phys_props, cell_chunk)
    for i in [1,3,5]
        @test phys_props.np[i,1] == 0.0
    end
    @test phys_props.np[4,1] == 1.0
    @test phys_props.n[4,1] == 4.0

    # even though we're not computing cell 2 properties, they are not reset either!
    @test phys_props.np[2,1] == 3.0
    @test phys_props.n[2,1] == 6.0

    # check that computation of number density works as well
    phys_props = PhysProps(pia; ndens_not_Np=true)

    cell_chunk = 1:3
    compute_props_sorted!(particles, pia, species_data, phys_props, grid, cell_chunk)
    for i in [1,3,4,5]
        @test phys_props.np[i,1] == 0.0
    end
    @test phys_props.np[2,1] == 3.0
    @test phys_props.n[2,1] == 6.0 / 2.0  # cell volume is 2.0


    cell_chunk = 4:4
    compute_props_sorted!(particles, pia, species_data, phys_props, grid, cell_chunk)
    for i in [1,3,5]
        @test phys_props.np[i,1] == 0.0
    end
    @test phys_props.np[2,1] == 3.0
    @test phys_props.n[2,1] == 6.0 / 2.0  # cell volume is 2.0
    @test phys_props.np[4,1] == 1.0
    @test phys_props.n[4,1] == 4.0 / 2.0  # cell volume is 2.0

    # 4 surf elements, 2 species
    surf_props1 = SurfProps(4, 2, [1.0, 1.0, 1.0, 1.0],
                            [1.0 0.0 0.0; -1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 -1.0 0.0;]')
    surf_props2 = SurfProps(4, 2, [1.0, 1.0, 1.0, 1.0],
                            [1.0 0.0 0.0; -1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 -1.0 0.0;]')
    surf_propst = SurfProps(4, 2, [1.0, 1.0, 1.0, 1.0],
                            [1.0 0.0 0.0; -1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 -1.0 0.0;]')

    surf_props1.flux_incident[:,1] .= [20.0, 10.0, 5.0, 17.0]
    surf_props1.flux_incident[:,2] .= [200.0, 100.0, 50.0, 170.0]

    surf_props2.flux_incident[:,1] .= [120.0, 110.0, 15.0, 117.0]
    surf_props2.flux_incident[:,2] .= [1200.0, 1100.0, 150.0, 1170.0]


    surf_props1.force[2,:,1] .= [71.0, 80.0, 99.0, 44.0]
    surf_props1.force[2,:,2] .= [-11.0, 19.0, 24.0, 33.0]

    surf_props2.force[2,:,2] .= [2.0, 4.0, 5.0, -6.0]
    surf_props2.force[3,:,1] .= [520.0, 190.0, 415.0, 17.0]
    surf_props2.force[3,:,2] .= [1705.0, 160.0, -152.0, 18.0]

    reduce_surf_props!(surf_propst, [surf_props1, surf_props2])
    @test maximum(abs.(surf_propst.flux_incident[:,1] - [140.0, 120.0, 20.0, 134.0])) < 2*eps()
    @test maximum(abs.(surf_propst.flux_incident[:,2] - [1400.0, 1200.0, 200.0, 1340.0])) < 2*eps()

    @test maximum(abs.(surf_propst.force[1,:,1] .- 0.0)) < 2*eps()
    @test maximum(abs.(surf_propst.force[1,:,2] .- 0.0)) < 2*eps()

    @test maximum(abs.(surf_propst.force[2,:,1] - [71.0, 80.0, 99.0, 44.0])) < 2*eps()
    @test maximum(abs.(surf_propst.force[2,:,2] - ([-11.0, 19.0, 24.0, 33.0] + [2.0, 4.0, 5.0, -6.0]))) < 2*eps()

    @test maximum(abs.(surf_propst.force[3,:,1] - [520.0, 190.0, 415.0, 17.0])) < 2*eps()
    @test maximum(abs.(surf_propst.force[3,:,2] - [1705.0, 160.0, -152.0, 18.0])) < 2*eps()

    # test sampling in a chunk
    particles = [ParticleVector(50)]

    pia = ParticleIndexerArray(grid.n_cells, 1)

    # sample only in cells 2, 3, and test that all other cells are empty
    # V = 2.0
    # 10 ppc: Np = 2e10 Fnum = 2e10 / 10 = 2e9
    cell_chunk = [2,3]
    sample_particles_equal_weight!(rng, grid, particles[1],
                                   pia, 1, species_data, 1e10, 300.0, 2e9, cell_chunk)
    compute_props_sorted!(particles, pia, species_data, phys_props, grid)

    @test phys_props.np[1,1] == 0
    @test phys_props.np[2,1] == 10
    @test phys_props.np[3,1] == 10
    @test phys_props.np[4,1] == 0
    @test phys_props.np[5,1] == 0
    @test phys_props.n[1,1] == 0
    @test phys_props.n[2,1] == 1e10
    @test phys_props.n[3,1] == 1e10
    @test phys_props.n[4,1] == 0
    @test phys_props.n[5,1] == 0

    # then sample in cell 4 and test that cells 2,3 unaffected
    # sample 20 particles with smaller weights
    cell_chunk = [4]
    sample_particles_equal_weight!(rng, grid, particles[1],
                                   pia, 1, species_data, 1e10, 300.0, 1e9, cell_chunk)
    compute_props_sorted!(particles, pia, species_data, phys_props, grid)

    @test phys_props.np[1,1] == 0
    @test phys_props.np[2,1] == 10
    @test phys_props.np[3,1] == 10
    @test phys_props.np[4,1] == 20
    @test phys_props.np[5,1] == 0
    @test phys_props.n[1,1] == 0
    @test phys_props.n[2,1] == 1e10
    @test phys_props.n[3,1] == 1e10
    @test phys_props.n[4,1] == 1e10
    @test phys_props.n[5,1] == 0
end