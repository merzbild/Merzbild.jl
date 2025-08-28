# Merzbild.jl

**Merzbild.jl** is a Direct Simulation Monte Carlo (DSMC) code written purely in Julia.
It provides all the necessary building blocks for building a DSMC simulation, i.e.
particle indexing, collisions, file I/O. Combining these blocks together is left up
to the user; examples can be found in the `simulations` directory.
It supports variable-weight DSMC simulations and ionized flow simulations,
PIC and Stochastic Fokker-Planck capabilities are in development.

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
For now, **Merzbild.jl** needs to be cloned to be run. Once cloned, navigate to the directory, run
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

More specific examples (ionizing collisions, specific merging algorithms) will be published later in a Tutorials section.

Finally, a full API reference is present, split into the [Merzbild.jl public API reference](@ref)
and [Merzbild.jl internal API reference](@ref).

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