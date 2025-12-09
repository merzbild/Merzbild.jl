@testset "merging_grid_merging" begin


    function create_particle_in_octant(octant, v_val; w=1.0)
        # assume octants symmetric around (0, 0, 0)
        # v_val > 0
        if octant >= 5
            v_z = v_val
        else
            v_z = -v_val
        end
        if octant % 2 == 1
            v_x = -v_val
        else
            v_x = v_val
        end
        if (octant == 3) || (octant == 4) || (octant == 7) || (octant == 8)
            v_y = v_val
        else
            v_y = -v_val
        end

        return Particle(w, [v_x, v_y, v_z], [1.0, -10.0, 3.0])
    end
    
    function create_24_3particles_in_octant(; weights=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    # create 3 particles in octant, each with weight == octant * weights[octant]
    # and velocity = 9.0 - octant - 0.5 / 9.0 - octant + 0.5 / 9.0 - octant
        vp = ParticleVector(24)

        i = 0
        for octant in 1:8
            i += 1
            Merzbild.update_particle_buffer_new_particle!(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant - 0.5, w=octant*weights[octant])
            i += 1
            Merzbild.update_particle_buffer_new_particle!(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant + 0.5, w=octant*weights[octant])
            i += 1
            Merzbild.update_particle_buffer_new_particle!(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant, w=octant*weights[octant])
        end

        return vp
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    rng = StableRNG(seed)

    Nx = 2
    Ny = 2
    Nz = 2

    phys_props::PhysProps = PhysProps(1, 1, [], Tref=1)

    Δabs = 2.5
    Δrel_xsmall = 5e-13

    n_particles = 5000
    particles::Vector{ParticleVector} = [ParticleVector(n_particles)]
    T0 = 300.0
    Fnum = 1e8
    vx0 = 2000.0
    vy0 = 500.0
    vz0 = -400.0
    pia = ParticleIndexerArray(0)

    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1,
                                   n_particles, species_data[1].mass, T0, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                   vx0=vx0, vy0=vy0, vz0=vz0)

    compute_props!(particles, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]

    mg = GridN2Merge(Nx, Ny, Nz, 3.5)

    merge_grid_based!(rng, mg, particles[1], pia, 1, 1, species_data, phys_props)

    @test pia.n_total[1] < n_particles
    @test pia.n_total[1] == pia.indexer[1,1].n_local

    compute_props!(particles, pia, species_data, phys_props)
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 3e-12
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 6.5e-12
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 3e-12
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 3e-12

    # test that empty cells are merged correctly
    weights_arr = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    particles24 = [create_24_3particles_in_octant(weights=weights_arr)]
    pia = ParticleIndexerArray(24)
    compute_props!(particles24, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]

    mg = GridN2Merge(2, 2, 2, 0.5)
    merge_grid_based!(rng, mg, particles24[1], pia, 1, 1, species_data, phys_props)
    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] < 24
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14

    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1

    @test e1 == pia.n_total[1]

    for i in s1:e1
        @test isnan(particles24[1][i].x[1]) == false
        @test particles24[1][i].w > 0
    end
    s2 = pia.indexer[1,1].start2
    e2 = pia.indexer[1,1].end2

    @test e2 == -1
    @test s2 == 0

    for i in s2:e2
        @test isnan(particles24[1][i].x[1]) == false
        @test particles24[1][i].w > 0
    end
    
    # now we split into 2 indexing groups
    weights_arr = [1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]
    particles24 = [create_24_3particles_in_octant(weights=weights_arr)]
    pia = ParticleIndexerArray(24)

    pia.indexer[1,1].n_group1 = 7
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 7

    pia.indexer[1,1].n_group2 = 17
    pia.indexer[1,1].start2 = 8
    pia.indexer[1,1].end2 = 24

    n_zero = 0

    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1
    s2 = pia.indexer[1,1].start2
    e2 = pia.indexer[1,1].end2

    for i in s1:e1
        if particles24[1][i].w == 0
            n_zero += 1
        end
    end
    for i in s2:e2
        if particles24[1][i].w == 0
            n_zero += 1
        end
    end
    @test n_zero == sum(weights_arr .== 0.0) * 3

    compute_props!(particles24, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]

    mg = GridN2Merge(2, 2, 2, 0.5)
    merge_grid_based!(rng, mg, particles24[1], pia, 1, 1, species_data, phys_props)
    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] < 24
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14

    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1
    s2 = pia.indexer[1,1].start2
    e2 = pia.indexer[1,1].end2

    @test e1 - s1 + 1 + e2 - s2 + 1 == pia.n_total[1]

    for i in s1:e1
        @test isnan(particles24[1][i].x[1]) == false
        @test particles24[1][i].w > 0
    end

    @test e2 > 0
    @test s2 > 0

    for i in s2:e2
        @test isnan(particles24[1][i].x[1]) == false
        @test particles24[1][i].w > 0
    end

    gm = GridN2Merge(4, [0.5, 0.7, 0.9])
    @test gm.Nx == 4
    @test gm.Ny == 4
    @test gm.Nz == 4
    @test maximum(abs.(gm.extent_multiplier .- [0.5, 0.7, 0.9])) < 2 * eps()

    gm = GridN2Merge(6, 1.5)
    @test gm.Nx == 6
    @test gm.Ny == 6
    @test gm.Nz == 6
    @test maximum(abs.(gm.extent_multiplier .- [1.5, 1.5, 1.5])) < 2 * eps()

    gm = GridN2Merge(4, 6, 2, 3.0, 2.0, 1.0)
    @test gm.Nx == 4
    @test gm.Ny == 6
    @test gm.Nz == 2
    @test maximum(abs.(gm.extent_multiplier .- [3.0, 2.0, 1.0])) < 2 * eps()
    # GridN2Merge(N::Int, extent_multiplier::T) where T <: AbstractArray = GridN2Merge(N, N, N, extent_multiplier)
    # GridN2Merge(N::Int, extent_multiplier::Float64) = GridN2Merge(N, N, N, extent_multiplier)
    # GridN2Merge(Nx::Int, Ny::Int, Nz::Int,
    #         extent_multiplier_x::Float64,
    #         extent_multiplier_y::Float64,
    #         extent_multiplier_z::Float64) = GridN2Merge(Nx, Ny, Nz, [extent_multiplier_x,
    #                                                                  extent_multiplier_y,
    #                                                                  extent_multiplier_z])
end