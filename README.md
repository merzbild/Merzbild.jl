# Merzbild.jl
**Merzbild.jl** is a work-in-progress DSMC code fully written in Julia,
designed to provide all the necessary components to build your own simulations.
This means that things like implementing a time loop over timesteps are left to the end user,
and the code provides only functionality such as
particle sampling and indexing, collisions, merging, computation of physical properties and I/O.

## Installation

Create a project where you want to use **Merzbild.jl**, run `julia --project=.`, type `]` to
enter the package manager, and run `add https://git-ce.rwth-aachen.de/georgii.oblapenko/merzbild.jl`.

## Usage
Currently, the way to use the code is to 1) clone it 2) create a new file in the `simulations` directory
3) add `include("../../../src/merzbild.jl")` and `using ..Merzbild` to the file.

Some usage examples can be found in the `simulations` directory.
A detailed overview of the structures required can be found in the documentation (upcoming),
but a short code snippet is given below. It simulates the relaxation of two gases initialized at different temperatures
to equilibrium.

```julia
include("../../../src/merzbild.jl")
using ..Merzbild

function run(seed)
    Random.seed!(seed)
    rng = Xoshiro(seed)

    # load species and interaction data
    species_list = load_species_list("path/to/data/particles.toml", ["Ar", "He"])
    interaction_data = load_interaction_data("data/vhs.toml", species_list)
    n_species = length(species_list)

    n_t = 800 # set number of timesteps

    # set numbers of particles/number densities and initial temperatures
    n_particles_Ar = 400
    n_particles_He = 4000
    T0_Ar = 3000.0
    T0_He = 360.0  # equilibrium T = 600.0K
    T0_list = [T0_Ar, T0_He]
    Fnum = 5e12

    # set timestep
    Δt = 2.5e-3

    # set simulation volume (we're doing a 0D simulation)
    V = 1.0

    # sample particles
    particles = [Vector{Particle}(undef, n_particles_Ar), Vector{Particle}(undef, n_particles_He)]
    sample_particles_equal_weight!(rng, particles[1], n_particles_Ar, T0_Ar, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    sample_particles_equal_weight!(rng, particles[2], n_particles_He, T0_He, species_list[2].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    # create struct for particle indexing
    particle_indexer = create_particle_indexer_array([n_particles_Ar, n_particles_He])

    # create struct for computation of physical properties
    phys_props = create_props(1, 2, [], Tref=T0_Ar)
    compute_props!(phys_props, particle_indexer, particles, species_list)

    # print properties at t=0
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)

    # create struct for output to netCDF file
    ds = create_netcdf_phys_props("output_multi_species.nc", phys_props, species_list)
    write_netcdf_phys_props(ds, phys_props, 0)

    # set up collision structs
    collision_factors = create_collision_factors(n_species)
    collision_data = create_collision_data()
    estimate_sigma_g_w_max!(collision_factors, interaction_data, species_list, T0_list, Fnum)

    # start time loop
    for ts in 1:n_t
        for s2 in 1:n_species
            # collide particles
            for s1 in s2:n_species
                if (s1 == s2)
                    ntc!(s1, 1, rng, collision_factors[s1,s1],
                         particle_indexer, collision_data, interaction_data[s1,s1], particles[s1],
                         Δt, V)
                else
                    ntc!(s1, s2, 1, rng, collision_factors[s1,s2], particle_indexer,
                         collision_data, interaction_data[s1,s2], particles[s1], particles[s2],
                         Δt, V)
                end
            end
        end

        # compute and write physical properties
        compute_props!(phys_props, particle_indexer, particles, species_list)
        write_netcdf_phys_props(ds, phys_props, ts)
    end

    # close output file
    close(ds)
end

run(1234)
```

### Usage notes

For now, bound checking is not turned off (via the `@inbounds` macro), so simulations may benefit from running with `--check-bounds=no`.
Running with `-O3` might also speed up things.

## Testing

The tests try to cover most of the functionality implemented in the code. They can be run by invoking `julia --project=. test/runtests.jl`.
Addition of test coverage (via `Coverage.jl`) is planned.