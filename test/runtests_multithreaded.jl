# Tests that actually use Threads.@threads. Unlike the serial-but-chunked tests in
# runtests_threading.jl, these only exercise anything when the process has >1 thread,
# so when the suite runs single-threaded each file is relaunched in a subprocess that
# has threads. Base.julia_cmd() propagates --check-bounds and --code-coverage.
const MULTITHREADED_TESTS = ["test_couette_chunked_threaded.jl"]
const MULTITHREADED_NTHREADS = min(4, Sys.CPU_THREADS)

if Threads.nthreads() > 1
    for testfile in MULTITHREADED_TESTS
        include(testfile)
    end
elseif MULTITHREADED_NTHREADS > 1
    for testfile in MULTITHREADED_TESTS
        @testset "$(testfile) (subprocess, $(MULTITHREADED_NTHREADS) threads)" begin
            cmd = `$(Base.julia_cmd()) --threads=$(MULTITHREADED_NTHREADS)
                   --project=$(dirname(Base.active_project()))
                   $(joinpath(@__DIR__, "run_threaded_testfile.jl")) $(testfile)`
            @test success(pipeline(cmd; stdout=stdout, stderr=stderr))
        end
    end
else
    @warn "Skipping multi-threaded tests: only one CPU thread available"
end
