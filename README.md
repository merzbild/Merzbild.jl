[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://merzbild.github.io/Merzbild.jl/dev/)
![Tests](https://github.com/merzbild/Merzbild.jl/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/merzbild/Merzbild.jl/branch/main/graph/badge.svg?token=MC2A2MM18K)](https://codecov.io/gh/merzbild/Merzbild.jl)
[![License: MPL2.0](https://img.shields.io/badge/License-MPL_2.0-success.svg)](https://opensource.org/license/mpl-2-0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14503197.svg)](https://doi.org/10.5281/zenodo.14503197)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Merzbild.jl
**Merzbild.jl** is a work-in-progress **variable-weight** DSMC code fully written in Julia,
designed to provide all the necessary components to build your own simulations.
This means that things like implementing a time loop over timesteps are left to the end user,
and the code provides only functionality such as
particle sampling and indexing, collisions, merging, computation of physical properties and I/O.
Currently the code supports serial and multithreaded operation.

## Installation
For now, **Merzbild.jl** needs to be cloned to be run. Once cloned, navigate to the directory, run `julia --project=.`, and in the
Julia interpreter run `using Pkg; Pkg.resolve(); Pkg.instantiate()` to install the required packages.
Running `Pkg.test()` afterwards will install the test environment dependencies and run the tests.
The package has been tested with the latest stable version of Julia 1.* (currently `1.11`), as well as the latest
Julia LTS version (currently `1.10`).

## Usage
Currently, the way to use the code is to 1) clone it 2) create a new file in the `simulations` directory
3) add `include("path/to/src/merzbild.jl")` and `using ..Merzbild` to the file.
A simulation file can be run by calling `julia --project=. simulations/path/to/simulation_file.jl` from the project
root directory.

Some usage examples can be found in the `simulations` directory.
A detailed overview of the structures required can be found in [the documentation](https://merzbild.github.io/Merzbild.jl/dev/),
but a short code snippet is given below. It simulates the relaxation of two gases initialized at different temperatures
to equilibrium.

```julia
# assuming the simulation file is directly in the simulations directory
include("../src/Merzbild.jl")
using ..Merzbild
using Random

function run(seed)
    Random.seed!(seed)
    rng = Xoshiro(seed)

    # load species and interaction data
    # path is correct if run from root directory of the Merzbild repo
    species_data = load_species_data("data/particles.toml", ["Ar", "He"])
    interaction_data = load_interaction_data("data/vhs.toml", species_data)
    n_species = length(species_data)

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

    # create struct for particle indexing
    pia = ParticleIndexerArray([0, 0])

    # sample particles
    particles::Vector{ParticleVector} = [ParticleVector(n_particles_Ar), ParticleVector(n_particles_He)]
    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, n_particles_Ar, species_data[1].mass, T0_Ar, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    sample_particles_equal_weight!(rng, particles[2], pia, 1, 2, n_particles_He, species_data[2].mass, T0_He, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    # create struct for computation of physical properties
    phys_props = PhysProps(pia)
    compute_props!(particles, pia, species_data, phys_props)

    # print properties at t=0
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)

    # create struct for output to netCDF file
    ds = NCDataHolder("output_multi_species.nc", species_data, phys_props)
    write_netcdf(ds, phys_props, 0)

    # set up collision structs
    collision_factors = create_collision_factors_array(n_species)
    collision_data = CollisionData()
    estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, Fnum)

    # start time loop
    for ts in 1:n_t
        for s2 in 1:n_species
            # collide particles
            for s1 in s2:n_species
                if (s1 == s2)
                    ntc!(rng, collision_factors[s1,s1,1], collision_data, interaction_data, particles[s1], pia, 1, s1, Δt, V)
                else
                    ntc!(rng, collision_factors[s1,s2,1], collision_data, interaction_data, particles[s1], particles[s2],
                         pia, 1, s1, s2, Δt, V)
                end
            end
        end

        # compute and write physical properties
        compute_props!(particles, pia, species_data, phys_props)
        write_netcdf(ds, phys_props, ts)
    end

    # print properties at last timestep
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)

    # close output file
    close_netcdf(ds)
end

run(1234)
```

### Usage notes
Simulations may benefit from running with `--check-bounds=no`, as some routines currently don't use `@inbounds`.
Running with `-O3` might also speed up things.

## Testing
The tests try to cover most of the functionality implemented in the code. They can be run by invoking by calling `using Pkg; Pkg.test()`.

## Speed
Detailed benchmarks can be found in [`BENCHMARKS.md`](BENCHMARKS.md).
For a serial 1-D Couette flow simulation, Merzbild.jl is up to 30% faster
than SPARTA.

## Citing
For now, the repository can be cited as
```bibtex
@misc{oblapenko2024merzbild,
  title={{M}erzbild.jl: A {J}ulia {DSMC} code},
  author={Oblapenko, Georgii},
  year={2024},
  month={12},
  howpublished={\url{https://github.com/merzbild/Merzbild.jl}},
  doi={10.5281/zenodo.14503197}
}
```

Depending on the specific functionality used, other citations may be warranted: please take a look
at the ["Overview of capabilities"](https://merzbild.github.io/Merzbild.jl/dev/overview_capabilities/) section in the documentation.

## Contributing
Please see [`CONTRIBUTING.MD`](CONTRIBUTING.MD) about some general guidelines on contributing to the development of Merzbild.jl.

## Acknowledgments
Dr. Georgii Oblapenko acknowledges the support of the German Research Foundation (DFG) via
the [SFB1481 research group](https://sfb1481.rwth-aachen.de).