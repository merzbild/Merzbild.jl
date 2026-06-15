using Aqua
using Test
using Merzbild
using Random
using NCDatasets
using SpecialFunctions
using StaticArrays
using StableRNGs
using LinearAlgebra
using ChunkSplitters

include("test_indexing.jl")  # first we test indexing routines
include("test_last_index_delete_add.jl")
include("test_constants.jl")  # test constants
include("test_species_data.jl")  # test loading of species data
include("test_computes.jl")  # we test functions that compute physical properties
include("test_io.jl")  # we test various I/O routines
include("test_sampling.jl")  # then we test sampling functions
include("test_grid_sampling.jl")  # then we test grid sampling functions
include("test_phase_box_sampling.jl")  # then we test phase box sampling
include("test_collision_utils.jl")  # then we test some collision utilities
include("test_collision_vhs.jl")  # then we test VHS collisions
include("test_collision_vhs_equal_weight.jl")  # then we test VHS collisions for equal-weight particles
include("test_electron_neutral_data.jl")  # test loading of XML LXCat-like data
include("test_acceleration.jl")  # test acceleration of charged particles
include("test_collision_ionizing.jl")  # test ionizing collisions
include("test_indexing_particlevector.jl")  # test particle vector and multi-cell indexing
include("test_indexing_particlebuffer.jl")  # test buffer in ParticleVector
include("test_pia_contiguous.jl")  # test squash_pia
include("test_debugging_functions.jl")  # test indexing debugging functions
include("test_grid_1D_uniform.jl")  # test sampling and computes on 1D uniform grid
include("test_grid_sorting.jl")  # test particle sorting routines
include("test_flux_computes.jl")  # test computation of fluxes in grid cells
include("test_io_fluxes.jl")  # test I/O of fluxes
include("test_convection_1D.jl")  # test particle convection and surface interaction on a 1D grid
include("test_collisions_1D.jl")  # test particle collisions on a grid
include("test_io_multidim.jl")  # test netCDF I/O of physical properties on a grid
include("test_io_particle.jl")  # test dump of particles
include("test_surface_props_1D_uniform.jl")  # test surface properties for 1-D uniform grid
include("test_io_surf.jl")  # test I/O of surface properties
include("test_misc.jl")  # various misc utility functions
include("test_collision_fp.jl")  # test Fokker-Planck collisions
include("test_1D_couette_fp.jl")  # test  1-D Couette flow, particle Fokker-Planck collisions
include("test_chunking.jl")  # test chunking, chunk-based computes, surface reduce, sampling on 1-D grid chunks
include("test_particle_exchange.jl")  # test swapping and pushing of particles in chunked simulations
include("test_particle_resort_after_exchange.jl")  # test re-sorting of particles after swap/push in chunked simulations
include("test_particle_index_sorting.jl")  # test re-sorting of particle indices
include("test_1_factorization.jl")  # test thread-safe 1-factorization of list of chunk pairs
include("test_collision_utils_swpm.jl")  # SWPM collision factors estimation
include("test_collisional_props_computes.jl")  # collisional property computes

include("runtests_merging.jl")  # tests for merging routines
include("runtests_reference_solutions.jl")  # tests that compare against reference solutions
include("runtests_malloc.jl")  # tests to check for unexpected memory allocations
# include("test_aqua.jl")


# tests assume that VHS data for Ar, He is
# ["Ar,Ar"]
# vhs_d = 4.11e-10
# vhs_o = 0.81
# vhs_Tref = 273.0
# comment = "Created via (1/2) mixing rule"

# ["Ar,He"]
# vhs_d = 3.25e-10
# vhs_o = 0.735
# vhs_Tref = 273.0
# comment = "Created via (1/2) mixing rule"

# ["He,He"]
# vhs_d = 2.33e-10
# vhs_o = 0.66
# vhs_Tref = 273.0
# comment = "Created via (1/2) mixing rule"