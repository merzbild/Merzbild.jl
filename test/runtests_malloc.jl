include("test_malloc_2species_equal_weight.jl")  # malloc with equal weight collisions and props computes
include("test_malloc_bkw_varweight_octree.jl")  # malloc with merging and variable-weight collisions
include("test_malloc_octree_1Dgrid.jl")  # malloc with octree merging on 1D grid
include("test_malloc_bkw_varweight_gridmerge.jl")  # malloc with grid merging and variable-weight collisions
include("test_malloc_gridmerge_1Dgrid.jl")  # malloc with grid merging on 1D grid
include("test_malloc_convection_1D.jl")  # malloc with convection on 1D grid
include("test_malloc_accelerate_ionize.jl")  # malloc with ionization reactions and acceleration


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