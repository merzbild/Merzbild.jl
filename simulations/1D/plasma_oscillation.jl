using Merzbild

"""
    run(L, nx, ppc, ndens, perturbation, n_timesteps, output_freq)

Run an electrostatic Particle-in-Cell simulation of a cold electron plasma oscillation in a
periodic 1-D domain. The electrons are initialized on a uniform lattice with a sinusoidal
perturbation of their positions on top of a fixed neutralizing ion background; they then oscillate
at the plasma frequency ``\\omega_p``. The computed physical properties are written to
`scratch/data/plasma_oscillation_[L]_[nx]_[ppc]_[ndens].nc`, the field energy is printed
to the standard output.

Positional arguments:
* `L`: length of domain in meters
* `nx`: number of grid cells
* `ppc`: number of electrons (and background ions) per cell
* `ndens`: electron number density
* `perturbation`: amplitude of the sinusoidal perturbation of the electron positions,
    relative to the domain length
* `n_timesteps`: total number of timesteps to run for
* `output_freq`: how frequently the output is written
"""
function run(L, nx, ppc, ndens, perturbation, n_timesteps, output_freq)
    # load particle data
    particles_data_path = joinpath(MERZBILD_DATA_PATH, "particles.toml")
    species_data = load_species_data(particles_data_path, ["e-", "He+"])

    # the timestep has to resolve the plasma oscillation: ω_p * Δt ≲ 0.2
    ω_p = plasma_frequency(1, species_data, ndens)
    Δt = 0.1 / ω_p
    println("ω_p = $ω_p, Δt = $Δt, ω_p * Δt = $(ω_p * Δt)")

    grid = Grid1DUniform(L, nx)

    # init particle vectors, particle indexer, grid particle sorter
    n_particles = ppc * nx
    particles = [ParticleVector(n_particles), ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 2)
    gridsorter = GridSortInPlace(grid, n_particles)

    # the electrons sit on a uniform lattice with a sinusoidal perturbation of their positions,
    # the ions are unperturbed and are never pushed, acting as a fixed neutralizing background
    w_particle = ndens * L / n_particles
    amplitude = perturbation * L
    k_wave = 2π / L

    for i in 1:n_particles
        x0 = L * (i - 0.5) / n_particles

        particles[1][i] = Particle(w_particle, [0.0, 0.0, 0.0], [x0 + amplitude * sin(k_wave * x0), 0.0, 0.0])
        particles[2][i] = Particle(w_particle, [0.0, 0.0, 0.0], [x0, 0.0, 0.0])
    end

    for species in 1:2
        pia.n_total[species] = n_particles
        pia.index_last[species] = n_particles
        pia.indexer[1,species].n_local = n_particles
        pia.indexer[1,species].n_group1 = n_particles
        pia.indexer[1,species].start1 = 1
        pia.indexer[1,species].end1 = n_particles
        pia.contiguous[species] = true

        sort_particles!(gridsorter, grid, particles[species], pia, species)
    end

    # the field quantities live on the nodes of the grid
    field_props = ElectrostaticFieldProps(grid)
    poisson_solver = PoissonSolver1DUniform(grid, PeriodicFieldBC1D(), PeriodicFieldBC1D())

    phys_props = PhysProps(pia)
    ds = NCDataHolder("scratch/data/plasma_oscillation_$(L)_$(nx)_$(ppc)_$(ndens).nc", species_data, phys_props)

    # leapfrog: the velocities are staggered by half a timestep with respect to the positions,
    # which is initialized by a single backward half-kick after the first field solve
    deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
    solve_poisson!(poisson_solver, field_props)
    accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, -0.5 * Δt)

    for t in 1:n_timesteps
        deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
        solve_poisson!(poisson_solver, field_props)

        accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, Δt)

        convect_particles_periodic!(grid, particles[1], pia, 1, Δt)
        sort_particles!(gridsorter, grid, particles[1], pia, 1)

        if t % output_freq == 0
            energy_field = 0.0
            for j in 1:nx
                energy_field += 0.5 * eps_0 * field_props.electric_field[j]^2 * grid.Δx
            end
            println("t = $t, ω_p * t = $(ω_p * t * Δt), field energy = $energy_field")

            compute_props_sorted!(particles, pia, species_data, phys_props)
            write_netcdf(ds, phys_props, t)
        end
    end

    close_netcdf(ds)
end

n_t = 2000
run(0.01, 32, 200, 1e14, 0.01, n_t, 10)
