include("test_merging_grid_indexing.jl")  # then we test basic indexing in grid-based merging
include("test_merging_grid_merging.jl")  # then we test grid-based merging in 0D
include("test_merging_grid_merging_2Dparticles.jl")  # then we test grid-based merging in 0D, 2D particles
include("test_merging_grid_merging_1Dparticles.jl")  # then we test grid-based merging in 0D, 1D particles
include("test_merging_grid_merging_0Dparticles.jl")  # then we test grid-based merging in 0D, 0D particles
include("test_merging_grid_buffer_sorting.jl")  # then we test that particles are freed correctly in NNLS-based merging
include("test_octree_bounds_and_splitting.jl")  # then we test bin bounds and splitting in octree merging
include("test_octree_sorting.jl")  # then we test sorting in octree merging
include("test_octree_merging.jl")  # then we test computation of props and octree merging
include("test_octree_merging_2Dparticles.jl")  # then we test computation of props and octree merging, Particle{2}
include("test_octree_merging_1Dparticles.jl")  # then we test computation of props and octree merging, Particle{1}
include("test_octree_merging_0Dparticles.jl")  # then we test computation of props and octree merging, Particle{0}
include("test_octree_merging_buffer_sorting.jl")  # then we test that particles are freed correctly in octree merging
include("test_nnls_utils.jl")  # test NNLS-based merging utils first
include("test_nnls_merging.jl")  # test NNLS-based merging
include("test_nnls_merging_buffer_sorting.jl")  # then we test that particles are freed correctly in NNLS-based merging
include("test_nnls_ratepreserving_merging.jl")  # test NNLS merging with approximate rate preservation
include("test_nnls_exact_ratepreserving_merging.jl")  # test NNLS merging with exact rate preservation 
include("test_roulette_merging.jl")  # roulette merge
include("test_merging_grid_merging_1D.jl")  # test grid-based merging in 1-D
include("test_octree_merging_1D.jl")  # test octree merging in 1-D

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