include("preamble.jl")  # shared imports, also used by the multi-threaded subprocess driver

include("runtests_basics.jl")  # fundamental blocks - indexing, sorting
include("runtests_computes.jl")  # various computes
include("test_acceleration.jl")  # test acceleration of charged particles
include("runtests_collisions_dsmc.jl")  # DSMC collisions
include("runtests_1DUniform.jl")  # 1D uniform grid, computes on grid, sorting routines
include("runtests_pic.jl")  # electrostatic Particle-in-Cell: Poisson solver, deposition, push
include("runtests_merging.jl")  # tests for merging routines
include("runtests_threading.jl")  # particle routines for multi-threading
include("runtests_multithreaded.jl")  # tests that actually run under Threads.@threads
include("runtests_FP.jl")  # tests for Fokker-Planck
include("runtests_io.jl")  # I/O tests
include("runtests_swpm.jl")  # SWPM tests
include("test_misc.jl")  # various misc utility functions
include("runtests_reference_solutions.jl")  # tests that compare against reference solutions
include("runtests_malloc.jl")  # tests to check for unexpected memory allocations
include("test_aqua.jl")


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