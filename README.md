# Merzbild.jl
**Merzbild.jl** is a work-in-progress DSMC code fully written in Julia, designed to provide all the necessary components to build your own simulations.
This means that things like implementing a time loop over timesteps are left to the end user, and the code provides only functionality such as
particle sampling and indexing, collisions, merging, computational of physical properties and I/O.

## Installation

## Usage

Some usage examples can be found in the `simulations` directory.

## Testing

The tests try to cover most of the functionality implemented in the code. They can be run by invoking `julia --project=. test/runtests.jl`.
Addition of test coverage (via `Coverage.jl`) is planned.