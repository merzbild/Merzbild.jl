include("../src/merging.jl")

using Test
using .Merging

include("test_computes.jl")  # first we test functions that compute physical properties
include("test_sampling.jl")  # then we test sampling functions
include("test_collision_utils.jl")  # then we test some collision utilities
include("test_bkw.jl")  # then we test the BKW relaxation case against reference solutions
include("test_2species.jl")  # then we test 2-species elastic collisions