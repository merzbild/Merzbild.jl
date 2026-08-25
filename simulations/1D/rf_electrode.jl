using Merzbild
using Random

"""
    run(seed, L, nx, ppc, ndens, T_e, T_i, V0, f_rf, n_timesteps, output_freq)

Run an electrostatic Particle-in-Cell simulation of a helium plasma between an RF-driven electrode:
a sine voltage applied to the left wall and a grounded wall on the right.
The particles are reflected specularly at both walls, so the particle number is conserved and no
surface charge accumulates. The computed physical properties are written to
`scratch/data/rf_electrode_[L]_[nx]_[ppc]_[V0]_[f_rf].nc`.

**Note**: this is currently unphysical, as no secondary emission occurs, electrons are not absorbed
at boundaries, and no collisions are modelled in the code.

Positional arguments:
* `seed`: random seed value
* `L`: length of domain in meters
* `nx`: number of grid cells
* `ppc`: number of particles per cell sampled per species
* `ndens`: number density of the electrons and of the ions
* `T_e`: initial electron temperature
* `T_i`: initial ion temperature
* `V0`: amplitude of the voltage applied to the left electrode
* `f_rf`: frequency of the voltage applied to the left electrode
* `n_timesteps`: total number of timesteps to run for
* `output_freq`: how frequently the output is written
"""
function run(seed, L, nx, ppc, ndens, T_e, T_i, V0, f_rf, n_timesteps, output_freq)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # load particle data
    particles_data_path = joinpath(MERZBILD_DATA_PATH, "particles.toml")
    species_data = load_species_data(particles_data_path, ["e-", "He+"])

    # the timestep has to resolve the electron plasma oscillation, the cell size the Debye length
    ω_p = plasma_frequency(1, species_data, ndens)
    λ_D = debye_length(ndens, T_e)
    Δt = 0.1 / ω_p

    grid = Grid1DUniform(L, nx)
    println("ω_p * Δt = $(ω_p * Δt), Δx / λ_D = $(grid.Δx / λ_D)")

    # init particle vectors, particle indexer, grid particle sorter
    n_particles = ppc * nx
    particles = [ParticleVector(n_particles), ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 2)
    gridsorter = GridSortInPlace(grid, n_particles)

    # sample particles
    # Fnum * ppc = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc
    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1, species_data, ndens, T_e, Fnum)
    sample_particles_equal_weight!(rng, grid, particles[2], pia, 2, species_data, ndens, T_i, Fnum)

    # particles are reflected specularly at both walls
    bc_list = (FullySpecularBC1D(), FullySpecularBC1D())

    # the left boundary is a driven electrode, the right one a grounded wall
    bc_left = DirichletFieldBC1D(0.0)
    bc_right = DirichletFieldBC1D(0.0)

    field_props = ElectrostaticFieldProps(grid)
    poisson_solver = PoissonSolver1DUniform(grid, bc_left, bc_right)

    phys_props = PhysProps(pia)
    ds = NCDataHolder("scratch/data/rf_electrode_$(L)_$(nx)_$(ppc)_$(V0)_$(f_rf).nc", species_data, phys_props)

    # leapfrog: initialize the half-step staggering of the velocities
    deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
    solve_poisson!(poisson_solver, field_props)

    for species in 1:2
        accelerate_electric_field_x!(grid, particles[species], pia, species, species_data, field_props, -0.5 * Δt)
    end

    for t in 1:n_timesteps
        # drive the electrode; only the right-hand side of the Poisson system depends on the value,
        # so no re-factorization of the matrix is needed
        bc_left.ϕ = V0 * sin(2π * f_rf * t * Δt)

        deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
        solve_poisson!(poisson_solver, field_props)

        for species in 1:2
            accelerate_electric_field_x!(grid, particles[species], pia, species, species_data, field_props, Δt)
            convect_particles!(rng, grid, bc_list, particles[species], pia, species, species_data, Δt)
            sort_particles!(gridsorter, grid, particles[species], pia, species)
        end

        if t % output_freq == 0
            println("t = $t, ϕ_electrode = $(bc_left.ϕ), max |E| = $(maximum(abs.(field_props.electric_field)))")

            compute_props_sorted!(particles, pia, species_data, phys_props)
            write_netcdf(ds, phys_props, t)
        end
    end

    close_netcdf(ds)
end

n_t = 2000
run(1234, 5e-3, 64, 200, 1e15, 11604.0, 300.0, 100.0, 13.56e6, n_t, 10)
