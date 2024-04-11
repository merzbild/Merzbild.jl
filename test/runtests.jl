include("../src/merzbild.jl")

using Test
using .Merzbild
using Random
using NCDatasets
using SpecialFunctions

include("test_computes.jl")  # first we test functions that compute physical properties
include("test_sampling.jl")  # then we test sampling functions
include("test_grid_sampling.jl")  # then we test grid sampling functions
include("test_collision_utils.jl")  # then we test some collision utilities
include("test_bkw.jl")  # then we test the BKW relaxation case against reference solutions
include("test_2species.jl")  # then we test 2-species elastic collisions
include("test_merging_grid_indexing.jl")  # then we test 2-species elastic collisions