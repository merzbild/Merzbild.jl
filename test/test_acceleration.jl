@testset "acceleration" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["He", "He+", "e-"])

    seed = 1234
    rng = StableRNG(seed)

    n_particles = 50000

    nv = 20

    n_dens_arr = [1e15, 1e10, 1e20]
    v0_arr = [[0.0, 15.0, 7.0],
              [-200.0, 35.0, 7.0],
              [1500.0, -2000.0, -75.5]]
    T0_arr = [400.0, 1000.0, 25000.0]

    index = 1
    n_sampled = [0, 0, 0]
    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles),
                                           Vector{Particle}(undef, n_particles),
                                           Vector{Particle}(undef, n_particles)]

    for (n_dens, v0, T0) in zip(n_dens_arr, v0_arr, T0_arr)
        n_sampled[index] = sample_maxwellian_on_grid!(rng, particles[index], nv, species_data[index].mass, T0, n_dens,
                                                      0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=3.5, cutoff_mult=8.0, noise=0.0, v_offset=v0)

        index += 1
    end

    pia = ParticleIndexerArray(n_sampled)
    pia2 = ParticleIndexerArray(n_sampled)

    n_offsets = [173, 188, 93]

    

    for i in 1:3
        @test n_offsets[i] < n_sampled[i]
        pia2.indexer[1,i].end1 = n_offsets[i]
        pia2.indexer[1,i].n_group1 = n_offsets[i]

        pia2.indexer[1,i].start2 = n_offsets[i] + 1
        pia2.indexer[1,i].end2 = n_sampled[i]
        pia2.indexer[1,i].n_group2 = n_sampled[i] - n_offsets[i]
    end

    phys_props::PhysProps = PhysProps(1, 3, [4, 6, 8], Tref=273.0)

    compute_props_with_total_moments!(particles, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,:]
    T0_computed = phys_props.T[1,:]
    vx_computed = phys_props.v[1,1,:]
    vy_computed = phys_props.v[2,1,:]
    vz_computed = phys_props.v[3,1,:]

    M4_computed = phys_props.moments[1,1,:]
    M6_computed = phys_props.moments[2,1,:]
    M8_computed = phys_props.moments[3,1,:]

    E = 997.0
    Δt = 1e-3

    for species in 1:3
        accelerate_constant_field_x!(particles[species], pia, 1, species, species_data, E, Δt)
    end

    # acceleration should not affect temperature/density/moments and y/z velocities
    # we will have some changes in T, M4/M6/M8 due to acceleration changing the x velocity
    # and subsequent change in round-off errors
    compute_props_with_total_moments!(particles, pia, species_data, phys_props)
    
    @test maximum(abs.(n0_computed - phys_props.n[1,:])) < eps()
    @test maximum(abs.(T0_computed - phys_props.T[1,:]) ./ T0_computed) < 1.5e-10
    @test maximum(abs.(vy_computed - phys_props.v[2,1,:])) < eps()
    @test maximum(abs.(vz_computed - phys_props.v[3,1,:])) < eps()
    @test maximum(abs.(M4_computed - phys_props.moments[1,1,:]) ./ M4_computed) < 2e-10
    @test maximum(abs.(M6_computed - phys_props.moments[2,1,:]) ./ M6_computed) < 3.3e-10
    @test maximum(abs.(M8_computed - phys_props.moments[3,1,:]) ./ M8_computed) < 4.9e-10

    # acceleration (via electric field) of neutral species should not affect it
    @test maximum(abs.(phys_props.v[1,1,1] - vx_computed[1])) < eps()

    acc_val = abs(E * Δt * Merzbild.q_e / 6.65e-27)
    @test abs.(phys_props.v[1,1,2] - (vx_computed[2] + E * Δt * Merzbild.q_e / 6.65e-27)) / acc_val < 1.5e-15

    acc_val = abs(E * Δt * Merzbild.q_e / 9.109383632e-31)
    @test maximum(abs.(phys_props.v[1,1,3] - (vx_computed[3] - E * Δt * Merzbild.q_e / 9.109383632e-31))) / acc_val < 1.5e-14

    # now we convect again, but this time the particle indexers have 2 groups of particles each
    for species in 1:3
        accelerate_constant_field_x!(particles[species], pia2, 1, species, species_data, E, Δt)
    end

    # acceleration should not affect temperature/density/moments and y/z velocities
    # we will have some changes in T, M4/M6/M8 due to acceleration changing the x velocity
    # and subsequent change in round-off errors
    compute_props_with_total_moments!(particles, pia, species_data, phys_props)
    
    @test maximum(abs.(n0_computed - phys_props.n[1,:])) < eps()
    @test maximum(abs.(T0_computed - phys_props.T[1,:]) ./ T0_computed) < 1.5e-10
    @test maximum(abs.(vy_computed - phys_props.v[2,1,:])) < eps()
    @test maximum(abs.(vz_computed - phys_props.v[3,1,:])) < eps()
    @test maximum(abs.(M4_computed - phys_props.moments[1,1,:]) ./ M4_computed) < 2e-10
    @test maximum(abs.(M6_computed - phys_props.moments[2,1,:]) ./ M6_computed) < 3.3e-10
    @test maximum(abs.(M8_computed - phys_props.moments[3,1,:]) ./ M8_computed) < 4.9e-10

    # acceleration (via electric field) of neutral species should not affect it
    @test maximum(abs.(phys_props.v[1,1,1] - vx_computed[1])) < eps()

    acc_val = abs(E * Δt * Merzbild.q_e / 6.65e-27)
    @test abs.(phys_props.v[1,1,2] - (vx_computed[2] + E * 2 * Δt * Merzbild.q_e / 6.65e-27)) / acc_val < 5e-14

    acc_val = abs(E * Δt * Merzbild.q_e / 9.109383632e-31)
    @test maximum(abs.(phys_props.v[1,1,3] - (vx_computed[3] - E * 2 * Δt * Merzbild.q_e / 9.109383632e-31))) / acc_val < 5e-14
end