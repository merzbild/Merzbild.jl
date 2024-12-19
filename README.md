[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://merzbild.github.io/Merzbild.jl/dev/)
![Tests](https://github.com/merzbild/Merzbild.jl/actions/workflows/ci.yml/badge.svg)
[![License: MPL2.0](https://img.shields.io/badge/License-MPL_2.0-success.svg)](https://opensource.org/license/mpl-2-0)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Merzbild.jl
**Merzbild.jl** is a work-in-progress DSMC code fully written in Julia,
designed to provide all the necessary components to build your own simulations.
This means that things like implementing a time loop over timesteps are left to the end user,
and the code provides only functionality such as
particle sampling and indexing, collisions, merging, computation of physical properties and I/O.

## Installation

For now, **Merzbild.jl** needs to be cloned to be run. Once cloned, navigate to the directory, run `julia --project=.`, and in the
Julia interpreter run `using Pkg; Pkg.resolve(); Pkg.instantiate()` to install the required packages.
Running `Pkg.test()` afterwards will install the test environment dependencies and run the tests.

## Usage
Currently, the way to use the code is to 1) clone it 2) create a new file in the `simulations` directory
3) add `include("path/to/src/merzbild.jl")` and `using ..Merzbild` to the file.

Some usage examples can be found in the `simulations` directory.
A detailed overview of the structures required can be found in the documentation (upcoming),
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
    write_netcdf_phys_props(ds, phys_props, 0)

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
        write_netcdf_phys_props(ds, phys_props, ts)
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

For now, bound checking is not turned off (via the `@inbounds` macro) except for the convection and particle sorting routines, so simulations may benefit from running with `--check-bounds=no`.
Running with `-O3` might also speed up things.

## Testing

The tests try to cover most of the functionality implemented in the code. They can be run by invoking by calling `using Pkg; Pkg.test()`.

## Speed
Comparing to SPARTA running in serial mode computing a Couette flow with 50000 particles and 50 cells (averaging over 36k timesteps after t>14000). Julia 1.9.4, with `--check-bounds=no -O3`.

Merzbild.jl:
```
 ──────────────────────────────────────────────────────────────────────────
                                  Time                    Allocations      
                         ───────────────────────   ────────────────────────
    Tot / % measured:         27.6s /  96.2%            699MiB /   0.4%    

 Section         ncalls     time    %tot     avg     alloc    %tot      avg
 ──────────────────────────────────────────────────────────────────────────
 sort             50.0k    8.90s   33.5%   178μs     0.00B    0.0%    0.00B
 convect          50.0k    6.70s   25.2%   134μs     0.00B    0.0%    0.00B
 collide          2.50M    5.80s   21.8%  2.32μs     0.00B    0.0%    0.00B
 props compute    36.0k    5.19s   19.5%   144μs     0.00B    0.0%    0.00B
 I/O                 51   4.44ms    0.0%  87.0μs   12.0KiB    0.4%     240B
 sampling             1   2.53ms    0.0%  2.53ms   3.05MiB   99.6%  3.05MiB
 ──────────────────────────────────────────────────────────────────────────
```

SPARTA:
```
Loop time of 31.3989 on 1 procs for 50000 steps with 50000 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 7.935      | 7.935      | 7.935      |   0.0 | 25.27
Coll    | 11.904     | 11.904     | 11.904     |   0.0 | 37.91
Sort    | 2.8201     | 2.8201     | 2.8201     |   0.0 |  8.98
Comm    | 0.0034597  | 0.0034597  | 0.0034597  |   0.0 |  0.01
Modify  | 8.7341     | 8.7341     | 8.7341     |   0.0 | 27.82
Output  | 0.00060415 | 0.00060415 | 0.00060415 |   0.0 |  0.00
Other   |            | 0.001508   |            |       |  0.00
```

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

Depending on the specific functionality used, other citations may be warranted (please look at the "Overview of capabilities" section
in the documentation).

## Contributing
Please see `CONTRIBUTING.MD` about some general guidelines on contributing to the development of Merzbild.jl

## Acknowledgments
Dr. Georgii Oblapenko acknowledges the support of the German Research Foundation (DFG) via
the [SFB1481 research group](https://sfb1481.rwth-aachen.de).