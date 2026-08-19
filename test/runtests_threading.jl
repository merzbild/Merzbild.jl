include("test_chunking.jl")  # test chunking, chunk-based computes, surface reduce, sampling on 1-D grid chunks
include("test_particle_exchange.jl")  # test swapping and pushing of particles in chunked simulations
include("test_occupancy_bounds.jl")  # test per-chunk occupied-cell bounds used to reject chunk pairs
include("test_particle_resort_after_exchange.jl")  # test re-sorting of particles after swap/push in chunked simulations
include("test_couette_varweight_octree_chunking.jl")  # test serial but chunked variable-weight simulation
include("test_1_factorization.jl")  # test thread-safe 1-factorization of list of chunk pairs
include("test_rebalance.jl")  # test rebalancing
include("test_couette_varweight_octree_chunk_and_rebalance.jl")  # test serial but chunked variable-weight simulation