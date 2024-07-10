include("../src/merzbild.jl")

using Test
using .Merzbild
using Random
using NCDatasets
using SpecialFunctions
using StaticArrays

include("test_indexing.jl")  # first we test indexing routines
include("test_computes.jl")  # we test functions that compute physical properties
include("test_sampling.jl")  # then we test sampling functions
include("test_grid_sampling.jl")  # then we test grid sampling functions
include("test_collision_utils.jl")  # then we test some collision utilities
include("test_bkw.jl")  # then we test the BKW relaxation case against reference solutions
include("test_2species.jl")  # then we test 2-species elastic collisions
include("test_merging_grid_indexing.jl")  # then we test basic indexing in grid-based merging
include("test_merging_grid_merging.jl")  # then we test grid-based merging in 0D
include("test_bkw_varweight_grid.jl")  # then we test variable-weight NTC collisions+grid-based merging in 0D (BKW relaxation)
include("test_octree_bounds_and_splitting.jl")  # then we test bin bounds and splitting in octree merging
include("test_octree_sorting.jl")  # then we test sorting in octree merging
include("test_octree_merging.jl")  # then we test computation of props and octree merging
include("test_bkw_varweight_octree.jl")  # then we test variable-weight NTC collisions+octree merging in 0D (BKW relaxation)
include("test_nnls_merging.jl")  # test NNLS-based merging
include("test_electron_neutral_data.jl")  # test loading of XML LXCat-like data
include("test_acceleration.jl")  # test acceleration of charge particles