include("test_bkw.jl")  # then we test the BKW relaxation case against reference solutions
include("test_bkw_equal_weight.jl")  # then we test the BKW relaxation case for equal-weight particles against reference solutions
include("test_2species.jl")  # then we test 2-species elastic collisions
include("test_2species_equal_weight.jl")  # then we test 2-species elastic collisions for equal-weight particles
# include("test_bkw_varweight_grid.jl")  # then we test variable-weight NTC collisions+grid-based merging in 0D (BKW relaxation)
# include("test_bkw_varweight_octree.jl")  # then we test variable-weight NTC collisions+octree merging in 0D (BKW relaxation)
# include("test_2species_varweight_octree.jl")  # 2-species variable-weight elastic collisions and merging
include("test_bkw_varweight_nnls.jl")  # test NNLS-based merging
# include("test_ionization_sim.jl")  # test 0D ionization simulation without event splitting
# include("test_ionization_sim_indexing.jl")  # test 0D ionization simulation with multiple neutral/ion species
include("test_1D_couette.jl")  # test Couette flow
# include("test_1D_couette_varweight.jl")  # test 1-D Couette flow, variable-weight DSMC, surface computes
include("test_1D_couette_nnls.jl")  # test 1-D Couette flow, variable-weight DSMC, NNLS merging surface computes
include("test_1D_couette_fp.jl")  # test  1-D Couette flow, particle Fokker-Planck collisions
# include("test_couette_varweight_octree_chunking.jl")  # test serial but chunked variable-weight simulation
# include("test_1D_couette_varweight_index_resort.jl")  # test 1-D Couette flow, variable-weight DSMC, index re-sorting
# include("test_bkw_varweight_octree_swpm.jl")  # SWPM for BKW test case
# include("test_1D_couette_varweight_swpm.jl")  # SWPM for Couette flow
include("test_2species_equal_weight0Dparticles.jl")  # test fixed-weight collisions with particles with dim(x)=
# include("test_bkw_varweight_octree_0Dparticle.jl")  # test variable weight BKW + octree merging with 0D particles"
# include("test_1D_couette_varweight_1DParticle.jl")  # test variable weight Couette + octree merging with 1D particles
# include("test_1D_couette_varweight_last_index.jl")  # test variable weight Couette with last_index and less squashing


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