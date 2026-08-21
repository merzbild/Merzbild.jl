include("test_collision_utils.jl")  # then we test some collision utilities
include("test_scattering_models.jl")  # elastic scattering model selection and per-pair dispatch
include("test_vhs_scattering.jl")  # then we test VHS collisions
include("test_vss_scattering.jl")  # then we test VSS and hard sphere scattering
include("test_collision_vhs.jl")  # then we test VHS collisions
include("test_collision_vhs_equal_weight.jl")  # then we test VHS collisions for equal-weight particles
include("test_electron_neutral_data.jl")  # test loading of XML LXCat-like data
include("test_collision_ionizing.jl")  # test ionizing collisions