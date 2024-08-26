@testset "collisions on 1D grid" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "pseudo_maxwell.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)
    np = 2000
    Fnum = 1e15
    T = 750.0

    grid = Grid1DUniform(10.0, 5)

    particles = [ParticleVector(np)]

    pia = ParticleIndexerArray(grid.n_cells, 1)
    

    cell = 2
    sample_particles_equal_weight!(rng, particles[1], pia, cell, 1,
                                   np, species_data[1].mass, T, Fnum,
                                   grid.cells[cell].xlo, grid.cells[cell].xhi,
                                   0.0, 1.0,
                                   0.0, 1.0;
                                   distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)

    @test pia.n_total[1] == np
    @test pia.indexer[cell, 1].start1 == 1
    @test pia.indexer[cell, 1].end1 == np

    phys_props = PhysProps(grid.n_cells, 1, [], Tref=1)
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.n[1, 1] == 0.0
    @test phys_props.n[2, 1] == np * Fnum
    @test phys_props.n[3, 1] == 0.0
    @test phys_props.n[4, 1] == 0.0
    @test phys_props.n[5, 1] == 0.0

    T_init = phys_props.T[2, 1]
    v_init = phys_props.v[:, 2, 1]


    collision_factors = create_collision_factors_array(1, grid.n_cells)
    collision_data = CollisionData()

    for cell in 1:grid.n_cells
        collision_factors[1, 1, cell].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T, Fnum)
    end

    Δt = 1e-3
    for cell in 1:grid.n_cells
        ntc!(rng, collision_factors[1, 1, cell],
             collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)

        @test collision_factors[1, 1, cell].n1 == phys_props.np[cell, 1]
        @test collision_factors[1, 1, cell].n2 == phys_props.np[cell, 1]

        if cell == 2
            @test collision_factors[1, 1, cell].n_coll > 0
        else
            collision_factors[1, 1, cell].n_coll == 0
        end
    end

    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.n[1, 1] == 0.0
    @test phys_props.n[2, 1] == np * Fnum
    @test phys_props.n[3, 1] == 0.0
    @test phys_props.n[4, 1] == 0.0
    @test phys_props.n[5, 1] == 0.0

    @test phys_props.np[2, 1] == np
    @test sum(phys_props.np[:, 1]) == np

    @test abs(phys_props.T[2, 1] - T_init) < 5e-13
    @test abs(phys_props.v[1, 2, 1] - v_init[1]) < 2e-14
    @test abs(phys_props.v[2, 2, 1] - v_init[2]) < 2e-14
    @test abs(phys_props.v[3, 2, 1] - v_init[3]) < 2e-14
end