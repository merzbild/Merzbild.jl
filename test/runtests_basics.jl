include("test_indexing.jl")  # first we test indexing routines
include("test_last_index_delete_add.jl")  # delete and add and impact on .last_index
include("test_constants.jl")  # test constants
include("test_species_data.jl")  # test loading of species data
include("test_io.jl")  # we test various I/O routines
include("test_sampling.jl")  # then we test sampling functions
include("test_grid_sampling.jl")  # then we test grid sampling functions
include("test_phase_box_sampling.jl")  # then we test phase box sampling
include("test_indexing_particlevector.jl")  # test particle vector and multi-cell indexing
include("test_indexing_particlebuffer.jl")  # test buffer in ParticleVector
include("test_pia_contiguous.jl")  # test squash_pia
include("test_debugging_functions.jl")  # test indexing debugging functions
include("test_particle_index_sorting.jl")  # test re-sorting of particle indices