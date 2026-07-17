@testset "linear FP, variable-weight particles" begin

    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "pseudo_maxwell.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data)
    
    n_particles = 100

    collision_data_fp = CollisionDataFP()

    particles = [ParticleVector(n_particles)]

    pia = ParticleIndexerArray(0)

    n_dens = 1e8
    T0 = 1000.0
    v0 = [-2500.0, 1000.0, 3000.0]

    Fnum = n_dens / n_particles

    # first: fixed-weight particle test 
    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, n_particles,
                                       species_data[1].mass, T0, Fnum,  0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                       vx0=v0[1], vy0=v0[2], vz0=v0[3])
    
    phys_props::PhysProps = PhysProps(1, 1)
    compute_props!(particles, pia, species_data, phys_props)

    n0_c = phys_props.n[1,1]
    v0_c = copy(phys_props.v[:,1,1])
    T0_c = phys_props.T[1,1]

    fp_linear!(rng, collision_data_fp, interaction_data[1,1], particles[1], pia, 1, 1, species_data, 1.0, 1.0)

    compute_props!(particles, pia, species_data, phys_props)

    @test abs(phys_props.n[1,1] - n0_c) / n0_c < 1e-15
    @test abs(maximum((phys_props.v[:,1,1] - v0_c) ./ v0_c)) < 1e-14
    @test abs(phys_props.T[1,1] - T0_c) / T0_c < 1e-15

    # we will use this in next step
    v_20_stored = copy(particles[1][20].v)
    v_60_stored = copy(particles[1][60].v)

    # now we will split particle indexing into group1, group2
    particles = [ParticleVector(n_particles)]

    rng2 = StableRNG(seed)
    pia = ParticleIndexerArray(0)
    sample_particles_equal_weight!(rng2, particles[1], pia, 1, 1, n_particles,
                                   species_data[1].mass, T0, Fnum,  0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                   vx0=v0[1], vy0=v0[2], vz0=v0[3])

    pia.indexer[1,1].n_group1 = 40
    pia.indexer[1,1].end1 = 40
    pia.indexer[1,1].n_group2 = 60
    pia.indexer[1,1].start2 = 41
    pia.indexer[1,1].end2 = 100

    fp_linear!(rng2, collision_data_fp, interaction_data[1,1], particles[1], pia, 1, 1, species_data, 1.0, 1.0)

    compute_props!(particles, pia, species_data, phys_props)

    @test abs(phys_props.n[1,1] - n0_c) / n0_c < 1e-15
    @test abs(maximum((phys_props.v[:,1,1] - v0_c) ./ v0_c)) < 1e-14
    @test abs(phys_props.T[1,1] - T0_c) / T0_c < 1e-15

    @test abs(maximum(particles[1][20].v - v_20_stored)) < 1e-10
    @test abs(maximum(particles[1][60].v - v_60_stored)) < 1e-10


    # now variable-weight
    particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(0)
    rng3 = StableRNG(seed)

    sample_particles_phase_box_weighted!(rng3, particles[1], pia, 1, 1, n_particles, species_data[1].mass, T0, n_dens,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=2.0, vx0=v0[1], vy0=v0[2], vz0=v0[3])

    compute_props!(particles, pia, species_data, phys_props)

    n0_c = phys_props.n[1,1]
    v0_c = copy(phys_props.v[:,1,1])
    T0_c = phys_props.T[1,1]

    fp_linear!(rng3, collision_data_fp, interaction_data[1,1], particles[1], pia, 1, 1, species_data, 1.0, 1.0)

    compute_props!(particles, pia, species_data, phys_props)

    @test abs(phys_props.n[1,1] - n0_c) / n0_c < 1e-15
    @test abs(maximum((phys_props.v[:,1,1] - v0_c) ./ v0_c)) < 1e-14
    @test abs(phys_props.T[1,1] - T0_c) / T0_c < 1e-14

    # we will use this in next step
    v_20_stored = copy(particles[1][20].v)
    v_60_stored = copy(particles[1][60].v)

    # now we will split particle indexing into group1, group2
    particles = [ParticleVector(n_particles)]

    rng4 = StableRNG(seed)
    pia = ParticleIndexerArray(0)

    sample_particles_phase_box_weighted!(rng4, particles[1], pia, 1, 1, n_particles, species_data[1].mass, T0, n_dens,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=2.0, vx0=v0[1], vy0=v0[2], vz0=v0[3])

    pia.indexer[1,1].n_group1 = 40
    pia.indexer[1,1].end1 = 40
    pia.indexer[1,1].n_group2 = 60
    pia.indexer[1,1].start2 = 41
    pia.indexer[1,1].end2 = 100

    fp_linear!(rng4, collision_data_fp, interaction_data[1,1], particles[1], pia, 1, 1, species_data, 1.0, 1.0)

    compute_props!(particles, pia, species_data, phys_props)

    @test abs(phys_props.n[1,1] - n0_c) / n0_c < 1e-15
    @test abs(maximum((phys_props.v[:,1,1] - v0_c) ./ v0_c)) < 1e-14
    @test abs(phys_props.T[1,1] - T0_c) / T0_c < 1e-14

    @test abs(maximum(particles[1][20].v - v_20_stored)) < 1e-10
    @test abs(maximum(particles[1][60].v - v_60_stored)) < 1e-10

    # now we will split particle indexing into group1, group2
    # and use particles with dim=2
    # can't check particle data since random sampling is different
    particles = [ParticleVector{2}(n_particles)]

    rng5 = StableRNG(seed)
    pia = ParticleIndexerArray(0)

    sample_particles_phase_box_weighted!(rng5, particles[1], pia, 1, 1, n_particles, species_data[1].mass, T0, n_dens,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=2.0, vx0=v0[1], vy0=v0[2], vz0=v0[3])

    pia.indexer[1,1].n_group1 = 40
    pia.indexer[1,1].end1 = 40
    pia.indexer[1,1].n_group2 = 60
    pia.indexer[1,1].start2 = 41
    pia.indexer[1,1].end2 = 100

    compute_props!(particles, pia, species_data, phys_props)

    n0_c = phys_props.n[1,1]
    v0_c = copy(phys_props.v[:,1,1])
    T0_c = phys_props.T[1,1]

    fp_linear!(rng5, collision_data_fp, interaction_data[1,1], particles[1], pia, 1, 1, species_data, 1.0, 1.0)

    compute_props!(particles, pia, species_data, phys_props)

    @test abs(phys_props.n[1,1] - n0_c) / n0_c < 1e-15
    @test abs(maximum((phys_props.v[:,1,1] - v0_c) ./ v0_c)) < 1e-14
    @test abs(phys_props.T[1,1] - T0_c) / T0_c < 1e-14


    # now we will split particle indexing into group1, group2
    # and use particles with dim=1
    # can't check particle data since random sampling is different
    particles = [ParticleVector{1}(n_particles)]

    rng6 = StableRNG(seed)
    pia = ParticleIndexerArray(0)

    sample_particles_phase_box_weighted!(rng6, particles[1], pia, 1, 1, n_particles, species_data[1].mass, T0, n_dens,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=2.0, vx0=v0[1], vy0=v0[2], vz0=v0[3])

    pia.indexer[1,1].n_group1 = 40
    pia.indexer[1,1].end1 = 40
    pia.indexer[1,1].n_group2 = 60
    pia.indexer[1,1].start2 = 41
    pia.indexer[1,1].end2 = 100

    compute_props!(particles, pia, species_data, phys_props)

    n0_c = phys_props.n[1,1]
    v0_c = copy(phys_props.v[:,1,1])
    T0_c = phys_props.T[1,1]

    fp_linear!(rng6, collision_data_fp, interaction_data[1,1], particles[1], pia, 1, 1, species_data, 1.0, 1.0)

    compute_props!(particles, pia, species_data, phys_props)

    @test abs(phys_props.n[1,1] - n0_c) / n0_c < 1e-15
    @test abs(maximum((phys_props.v[:,1,1] - v0_c) ./ v0_c)) < 1e-14
    @test abs(phys_props.T[1,1] - T0_c) / T0_c < 1e-14
end