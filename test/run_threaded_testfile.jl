include(joinpath(@__DIR__, "preamble.jl"))

# driver used by runtests_multithreaded.jl to run a single test file in a subprocess
# started with more than one thread; a failing top-level @testset throws, giving a
# non-zero exit code that the parent process picks up
include(joinpath(@__DIR__, ARGS[1]))
