# Inspired by
# - https://gitlab.com/tkpapp/GitlabJuliaDemo.jl by Tamas K. Papp
# - https://github.com/tpapp/LocalCoverage.jl by Tamas K. Papp
# - https://github.com/trixi-framework/Trixi.jl/blob/main/test/coverage/coverage.jl by the Trixi.jl team

using Coverage

const report_dir = "coverage_report"
const lcov_info_file = "lcov.info"

# Change path to root directory
cd(joinpath(@__DIR__, "..", "..")) do
    # Process coverage files
    processed = process_folder("src")

    # Uncomment the following line once Codecov support is enabled
    # Codecov.submit_local(processed)

    # Calculate coverage
    covered_lines, total_lines = get_summary(processed)
    percentage = covered_lines / total_lines * 100

    # Print coverage in a format that can be easily parsed
    println("($(percentage)%) covered")

    # Try to generate a coverage report
    isdir(report_dir) || mkdir(report_dir)
    tracefile = joinpath(report_dir, lcov_info_file)
    Coverage.LCOV.writefile(tracefile, processed)

    # Clean up .cov files
    clean_folder("src")
    clean_folder("test")
end