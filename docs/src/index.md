# Merzbild.jl

**Merzbild.jl** is a Direct Simulation Monte Carlo (DSMC) code written purely in Julia.
It provides all the necessary building blocks for building a DSMC simulation, i.e.
particle indexing, collisions, file I/O. Combining these blocks together is left up
to the user; examples can be found in the `simulations` directory.
It supports variable-weight DSMC and stochastic linear Fokker-Planck simulations and ionized flow simulations;
PIC capabilities are in development.

The goals are to provide a modular, thoroughly tested, easy-to-read and easy-to-extend code
for quick implementation and testing of new ideas.

Serial and multithreaded simulations are possible.

## Brief overview of capabilities
Currently the code supports spatially homogeneous (0D) and 1D uniform grid simulations.
The table below lists support for the fixed- and variable-weight versions of the code.

|                        | **0D**                                        | **1D** |
|:----------------------:|:-----------------------------------------:|:----:|
| Fixed-weight DSMC      | ✅                                        | ✅ |
| Variable-weight DSMC   | ✅ | ✅ |
| Fixed-weight Fokker-Planck| Linear | Linear |

A more detailed overview of the capabilities is given on the [Overview of capabilities](@ref) page.

The output format is NetCDF4.

## Installation

### Package installation
You can either install **Merzbild.jl** as a package via the Julia package manager (`using Pkg; Pkg.add("Merzbild")`).
You can run tests by entering the package manager (`]`) and running `test Merzbild`.
An example simulation from the `simulations` folder can be run by typing
`include(joinpath(MERZBILD_SIMULATIONS_PATH, "relative/path/to_simulation.jl"))`.
Simulations bundled with Merzbild output data to `scratch/data` and will crash if the directory is not present.

The constant `MERZBILD_DATA_PATH` is exported and points to the `data` directory of the package which stores
species and interaction data.

### Development installation
To develop functionality within the package, clone the package. Once cloned, navigate to the directory, run
```
julia --project=.
```
and in the Julia interpreter run
```julia
using Pkg; Pkg.resolve(); Pkg.instantiate();
```
to install the required packages.
Running `Pkg.test()` afterwards will install the test environment dependencies and run the tests.

## Usage
Currently, the way to use the code is to
  1. clone it
  2. create a new file in the `simulations` directory
  3. add `include("path/to/src/merzbild.jl")` and `using ..Merzbild` to the file.

## Documentation
The documentation assumes a certain level of pre-existing knowledge of the DSMC approach.
Basic building blocks and operations (particle indexing, sampling, collisions, I/O) are covered in the
Getting Started section.

More specific examples and some advanced code aspects are found in the Tutorials section.
A list of various algorithms implemented in the code is given in the Implemented algorithms section.

Finally, a full API reference is present, split into the [Merzbild.jl public API reference](@ref)
and [Merzbild.jl internal API reference](@ref).

## Example simulations
Various example simulations are available in the `simulations/` directory; the path to the directory as bundled with the code is
exported as `MERZBILD_SIMULATIONS_PATH`.

## Citing
You can for now cite the repository as
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

Depending on the functionality used, other citations may be warranted, please look at the
[Overview of capabilities](@ref) page to see which algorithms and models have been implemented in Merzbild.