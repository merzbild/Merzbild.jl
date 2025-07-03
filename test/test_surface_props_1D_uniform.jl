@testset "1D uniform grid surface computations" begin
    Δlarge = 24.0

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    grid = Grid1DUniform(4.0, 8)

    ppc = 4

    particles = [ParticleVector(ppc * grid.n_cells)]

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    pia = ParticleIndexerArray(grid.n_cells, 1)

    surf_props = SurfProps(pia, grid)

    @test surf_props.n_elements == 2
    @test surf_props.n_species == 1
    @test surf_props.areas == [1.0, 1.0]
    @test surf_props.inv_areas == [1.0, 1.0]
    @test surf_props.normals[:,1] == [1.0, 0.0, 0.0]
    @test surf_props.normals[:,2] == [-1.0, 0.0, 0.0]
    @test size(surf_props.np) == (2,1)
    @test size(surf_props.flux_incident) == (2,1)
    @test size(surf_props.flux_reflected) == (2,1)
    @test size(surf_props.force) == (3,2,1)
    @test size(surf_props.normal_pressure) == (2,1)
    @test size(surf_props.shear_pressure) == (3,2,1)
    @test size(surf_props.kinetic_energy_flux) == (2,1)

    # 4 particles in cell 1, with velocity -i^2
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle(particles[1], i)

        w = i + 0.5
        particles[1][i] = Particle(w, [-i^2, i^2, 0], [0.25, 0.0, 0.0])
    end

    # 4 particles in cell 8, with velocity 10.0
    for i in 5:8
        Merzbild.update_particle_buffer_new_particle(particles[1], i)

        w = i + 2.5
        particles[1][i] = Particle(w, [10.0, 0, 3.0], [3.9, 0.0, 0.0])
    end

    pia.n_total[1] = 8
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4
    pia.indexer[1,1].n_group1 = 4

    pia.indexer[1,1].start2 = -1
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0

    pia.indexer[8,1].n_local = 4
    pia.indexer[8,1].start1 = 5
    pia.indexer[8,1].end1 = 8
    pia.indexer[8,1].n_group1 = 4

    pia.indexer[8,1].start2 = -1
    pia.indexer[8,1].end2 = -1
    pia.indexer[8,1].n_group2 = 0


    clear_props!(surf_props)

    # let's assume particle 3 and particle 8 have hit the L and R walls
    # particle, species, surf_props, surface_element_id, mass, Δt
    Merzbild.update_surface_incident!(particles[1][3], 1, surf_props, 1)
    Merzbild.update_surface_incident!(particles[1][8], 1, surf_props, 2)

    @test surf_props.np == [1.0;1.0;;]

    @test surf_props.flux_incident == [3.5;10.5;;]
    @test surf_props.flux_reflected == [0.0;0.0;;]

    @test abs(surf_props.force[1,1,1] - (3.5 * -9)) < 2*eps() # velocity is -9
    @test abs(surf_props.force[1,2,1] - (10.5 * 10)) < 2*eps() # velocity is 10

    @test abs(surf_props.force[2,1,1] - (3.5 * 9)) < 2*eps()  # y-velocity is 9
    @test surf_props.force[2,2,1] == 0.0
    
    @test surf_props.force[3,1,1] == 0.0
    @test abs(surf_props.force[3,2,1] - (10.5 * 3)) < 2*eps()  # z-velocity is 3

    p_normal_left = -9
    p_normal_right = -10

    @test abs(surf_props.normal_pressure[1,1] - (-p_normal_left * 3.5)) < 2*eps()
    @test abs(surf_props.normal_pressure[2,1] - (-p_normal_right * 10.5)) < 2*eps()

    @test abs(surf_props.shear_pressure[1,1,1] - (0 * 3.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[2,1,1] - (9 * 3.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[3,1,1] - (0 * 3.5)) < 2*eps()

    @test abs(surf_props.shear_pressure[1,2,1] - (0 * 10.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[2,2,1] - (0 * 10.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[3,2,1] - (3 * 10.5)) < 2*eps()

    @test abs(surf_props.kinetic_energy_flux[1,1] - (0.5 * 3.5 * (9^2 + 9^2))) < 2*eps()
    @test abs(surf_props.kinetic_energy_flux[2,1] - (0.5 * 10.5 * (10^2 + 3^2))) < 2*eps()

    # specular reflection of left particle, was [-9.0, 9.0, 0.0]
    particles[1][3].v = SVector{3,Float64}([9.0, 9.0, 0.0])

    # right particle reflected with completely new velocity, [-1.0, -2.0, 4.0]
    # but also reduced weight (assume part is absorbed or something)
    particles[1][8].w = 10.0
    particles[1][8].v = SVector{3,Float64}([-1.0, -2.0, -4.0])

    Merzbild.update_surface_reflected!(particles[1][3], 1, surf_props, 1)
    Merzbild.update_surface_reflected!(particles[1][8], 1, surf_props, 2)


    @test surf_props.np == [1.0;1.0;;]

    @test surf_props.flux_incident == [3.5;10.5;;]
    @test surf_props.flux_reflected == [-3.5;-10.0;;]
    
    # new forces
    @test abs(surf_props.force[1,1,1] - (3.5 * -9 - 3.5*9)) < 2*eps() # -9 -> 9
    @test abs(surf_props.force[1,2,1] - (10.5 * 10 - 10.0 * -1)) < 2*eps() # 10 -> -1

    @test abs(surf_props.force[2,1,1] - 0.0) < 2*eps()  # 9 -> 9
    @test abs(surf_props.force[2,2,1] - (-10.0 * -2)) < 2*eps()  # 0 -> -2
    
    @test surf_props.force[3,1,1] == 0.0  # 0 -> 0
    @test abs(surf_props.force[3,2,1] - (10.5 * 3 - 10.0 * -4)) < 2*eps()  # 3 -> -4

    p_normal_left_reflected = 9.0
    p_normal_right_reflected = 1.0

    @test abs(surf_props.normal_pressure[1,1] - (-p_normal_left * 3.5 + p_normal_left_reflected * 3.5)) < 2*eps()
    @test abs(surf_props.normal_pressure[2,1] - (-p_normal_right * 10.5 + p_normal_right_reflected * 10.0)) < 2*eps()

    @test abs(surf_props.shear_pressure[1,1,1] - (0 * 3.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[2,1,1] - (9 * 3.5 - 9 * 3.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[3,1,1] - (0 * 3.5)) < 2*eps()

    @test abs(surf_props.shear_pressure[1,2,1] - (0 * 10.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[2,2,1] - (0 * 10.5 - -2 * 10.0)) < 2*eps()
    @test abs(surf_props.shear_pressure[3,2,1] - (3 * 10.5 - -4 * 10.0)) < 2*eps()

    @test abs(surf_props.kinetic_energy_flux[1,1]) < 2*eps()
    @test abs(surf_props.kinetic_energy_flux[2,1] - (0.5 * 10.5 * (10^2 + 3^2) - 0.5 * 10.0 * (1^2 + 2^2 + 4^2))) < 2*eps()

    # make something != 1 so that we can test the actual scaling
    surf_props.areas = [2.0, 4.0]
    surf_props.inv_areas = [0.5, 0.25]
    Δt = 1e-20  # we make this super small so that any errors remain larger than machine precision
    Merzbild.surface_props_scale!(1, surf_props, species_data, Δt)

    f = species_data[1].mass * surf_props.inv_areas / Δt

    @test abs(surf_props.flux_incident[1,1] - f[1] * 3.5) < 2*eps()
    @test abs(surf_props.flux_incident[2,1] - f[2] * 10.5) < 2*eps()
    @test abs(surf_props.flux_reflected[1,1] - -f[1] * 3.5) < 2*eps()
    @test abs(surf_props.flux_reflected[2,1] - -f[2] * 10.0) < 2*eps()

    @test abs(surf_props.normal_pressure[1,1] - f[1]*(-p_normal_left * 3.5 + p_normal_left_reflected * 3.5)) < 2*eps()
    @test abs(surf_props.normal_pressure[2,1] - f[2]*(-p_normal_right * 10.5 + p_normal_right_reflected * 10.0)) < 2*eps()

    @test abs(surf_props.shear_pressure[1,1,1] - f[1]*(0 * 3.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[2,1,1] - f[1]*(9 * 3.5 - 9 * 3.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[3,1,1] - f[1]*(0 * 3.5)) < 2*eps()

    @test abs(surf_props.shear_pressure[1,2,1] - f[2]*(0 * 10.5)) < 2*eps()
    @test abs(surf_props.shear_pressure[2,2,1] - f[2]*(0 * 10.5 - -2 * 10.0)) < 2*eps()
    @test abs(surf_props.shear_pressure[3,2,1] - f[2]*(3 * 10.5 - -4 * 10.0)) < 2*eps()

    @test abs(surf_props.kinetic_energy_flux[1,1] * f[1]) < 2*eps()
    @test abs(surf_props.kinetic_energy_flux[2,1] - f[2]*(0.5 * 10.5 * (10^2 + 3^2) - 0.5 * 10.0 * (1^2 + 2^2 + 4^2))) < 2*eps()

    # reset back
    surf_props.areas = [1.0, 1.0]
    surf_props.inv_areas = [1.0, 1.0]

    # test clearing of data
    clear_props!(surf_props)
    
    # we don't over-write the fixed stuff
    @test surf_props.n_elements == 2
    @test surf_props.n_species == 1
    @test surf_props.areas == [1.0, 1.0]
    @test surf_props.inv_areas == [1.0, 1.0]
    @test surf_props.normals[:,1] == [1.0, 0.0, 0.0]
    @test surf_props.normals[:,2] == [-1.0, 0.0, 0.0]

    # the other things are set to 0
    @test maximum(abs.(surf_props.np)) == 0.0
    @test maximum(abs.(surf_props.flux_incident)) == 0.0
    @test maximum(abs.(surf_props.flux_reflected)) == 0.0
    @test maximum(abs.(surf_props.force)) == 0.0
    @test maximum(abs.(surf_props.normal_pressure)) == 0.0
    @test maximum(abs.(surf_props.shear_pressure)) == 0.0
    @test maximum(abs.(surf_props.kinetic_energy_flux)) == 0.0

    # test averaging of props
    # fill with dummy data
    surf_props.np[1,1] = 20
    surf_props.np[2,1] = 10

    surf_props.flux_incident[1,1] = 6
    surf_props.flux_reflected[2,1] = 11

    surf_props.force[1,1,1] = 1
    surf_props.force[2,1,1] = 3
    surf_props.force[3,1,1] = -8

    surf_props.force[1,2,1] = 4
    surf_props.force[2,2,1] = 9
    surf_props.force[3,2,1] = 5

    surf_props.normal_pressure[1,1] = 11
    surf_props.normal_pressure[2,1] = 9

    surf_props.shear_pressure[1,1,1] = 4
    surf_props.shear_pressure[2,1,1] = 5
    surf_props.shear_pressure[3,1,1] = 6

    surf_props.shear_pressure[1,2,1] = 11
    surf_props.shear_pressure[2,2,1] = 9
    surf_props.shear_pressure[3,2,1] = 2
    
    surf_props.kinetic_energy_flux[1,1] = 3
    surf_props.kinetic_energy_flux[2,1] = 7

    surf_props_avg = SurfProps(pia, grid)

    avg_props!(surf_props_avg, surf_props, 2)

    surf_props.np[1,1] = 10
    surf_props.np[2,1] = 20

    surf_props.flux_incident[1,1] = 5
    surf_props.flux_reflected[2,1] = 3

    surf_props.force[1,1,1] = 2
    surf_props.force[2,1,1] = 4
    surf_props.force[3,1,1] = 8

    surf_props.force[1,2,1] = 3
    surf_props.force[2,2,1] = 1
    surf_props.force[3,2,1] = 9

    surf_props.normal_pressure[1,1] = 12
    surf_props.normal_pressure[2,1] = 6

    surf_props.shear_pressure[1,1,1] = 1
    surf_props.shear_pressure[2,1,1] = 9
    surf_props.shear_pressure[3,1,1] = 8

    surf_props.shear_pressure[1,2,1] = 3
    surf_props.shear_pressure[2,2,1] = 4
    surf_props.shear_pressure[3,2,1] = 8
    
    surf_props.kinetic_energy_flux[1,1] = 2
    surf_props.kinetic_energy_flux[2,1] = 6

    avg_props!(surf_props_avg, surf_props, 2)

    @test abs(surf_props_avg.np[1,1] - 15) < 2 * eps()
    @test abs(surf_props_avg.np[2,1] - 15) < 2 * eps()

    @test abs(surf_props_avg.flux_incident[1,1] - 5.5) < 2 * eps()
    @test abs(surf_props_avg.flux_reflected[2,1] - 7) < 2 * eps()
 
    @test abs(surf_props_avg.force[1,1,1] - 1.5) < 2 * eps()
    @test abs(surf_props_avg.force[2,1,1] - 3.5) < 2 * eps()
    @test abs(surf_props_avg.force[3,1,1] - 0.0) < 2 * eps()

    @test abs(surf_props_avg.force[1,2,1] - 3.5) < 2 * eps()
    @test abs(surf_props_avg.force[2,2,1] - 5) < 2 * eps()
    @test abs(surf_props_avg.force[3,2,1] - 7) < 2 * eps()

    @test abs(surf_props_avg.normal_pressure[1,1] - 11.5) < 2 * eps()
    @test abs(surf_props_avg.normal_pressure[2,1] - 7.5) < 2 * eps()

    @test abs(surf_props_avg.shear_pressure[1,1,1] - 2.5) < 2 * eps()
    @test abs(surf_props_avg.shear_pressure[2,1,1] - 7) < 2 * eps()
    @test abs(surf_props_avg.shear_pressure[3,1,1] - 7) < 2 * eps()

    @test abs(surf_props_avg.shear_pressure[1,2,1] - 7) < 2 * eps()
    @test abs(surf_props_avg.shear_pressure[2,2,1] - 6.5) < 2 * eps()
    @test abs(surf_props_avg.shear_pressure[3,2,1] - 5) < 2 * eps()

    @test abs(surf_props_avg.kinetic_energy_flux[1,1] - 2.5) < 2 * eps()
    @test abs(surf_props_avg.kinetic_energy_flux[2,1] - 6.5) < 2 * eps()

    # test convection, contiguous
    # specular walls
    boundaries = MaxwellWalls1D(species_data, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)

    ppc = 4
    particles = [ParticleVector(ppc * grid.n_cells)]

    # 4 particles in cell 1, all of them reach the wall
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle(particles[1], i)

        w = 1.0
        particles[1][i] = Particle(w, [-1.0, 0, 0], [0.1, 0.0, 0.0])
    end

    # 4 particles in cell 8, only one of them reaches the wall
    for i in 5:8
        Merzbild.update_particle_buffer_new_particle(particles[1], i)

        w = 3.0
        vx = 0.0
        if i == 5
            vx = 1.0
        end
        particles[1][i] = Particle(w, [vx, 0, 3.0], [3.9, 0.0, 0.0])
    end

    pia.n_total[1] = 8
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4
    pia.indexer[1,1].n_group1 = 4

    pia.indexer[1,1].start2 = -1
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0

    pia.indexer[8,1].n_local = 4
    pia.indexer[8,1].start1 = 5
    pia.indexer[8,1].end1 = 8
    pia.indexer[8,1].n_group1 = 4

    pia.indexer[8,1].start2 = -1
    pia.indexer[8,1].end2 = -1
    pia.indexer[8,1].n_group2 = 0

    clear_props!(surf_props)

    convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, surf_props, 0.2)

    @test surf_props.np == [4.0;1.0;;]


    f = species_data[1].mass * surf_props.inv_areas / 0.2

    @test abs(surf_props.flux_incident[1,1]/f[1] - 4.0) < 3*eps()
    @test abs(surf_props.flux_incident[2,1]/f[2] - 3.0) < 3*eps()
    @test abs(surf_props.flux_reflected[1,1]/f[1] - -4.0) < 3*eps()
    @test abs(surf_props.flux_reflected[2,1]/f[2] - -3.0) < 3*eps()
    # we test the computation of the quantities before, so this should be sufficient to verify that 


    clear_props!(surf_props)
    # test convection, non-contiguous
    particles = [ParticleVector(8)]

    # 4 particles in cell 1, all of them reach the wall
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle(particles[1], i)

        w = 1.0
        particles[1][i] = Particle(w, [-1.0, 0, 0], [0.1, 0.0, 0.0])
    end

    # 4 particles in cell 8, only one of them reaches the wall
    for i in 5:8
        Merzbild.update_particle_buffer_new_particle(particles[1], i)

        w = 3.0
        vx = 0.0
        if i == 5
            vx = 1.0
        end
        particles[1][i] = Particle(w, [vx, 0, 3.0], [3.9, 0.0, 0.0])
    end
    
    
    pia.n_total[1] = 8
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4
    pia.indexer[1,1].n_group1 = 4

    pia.indexer[1,1].start2 = -1
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0

    pia.indexer[8,1].n_local = 4
    pia.indexer[8,1].start1 = 5
    pia.indexer[8,1].end1 = 8
    pia.indexer[8,1].n_group1 = 4

    pia.indexer[8,1].start2 = -1
    pia.indexer[8,1].end2 = -1
    pia.indexer[8,1].n_group2 = 0

    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    Merzbild.delete_particle_end!(particles[1], pia, 1, 1)
    pia.contiguous[1] = false

    convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, surf_props, 0.2)

    @test surf_props.np == [2.0;1.0;;]

    f = species_data[1].mass * surf_props.inv_areas / 0.2

    @test abs(surf_props.flux_incident[1,1]/f[1] - 2.0) < 3*eps()
    @test abs(surf_props.flux_incident[2,1]/f[2] - 3.0) < 3*eps()
    @test abs(surf_props.flux_reflected[1,1]/f[1] - -2.0) < 3*eps()
    @test abs(surf_props.flux_reflected[2,1]/f[2] - -3.0) < 3*eps()
end 