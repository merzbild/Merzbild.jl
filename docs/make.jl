using Documenter
using Merzbild

makedocs(
    sitename = "Merzbild",
    format = Documenter.HTML(; size_threshold=400_000),
    modules = [Merzbild],
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Overview of capabilities" => "overview_capabilities.md",
        "Getting started" => [
            "Overview of basic building blocks" => "overview_blocks.md",
            "Fixed-weight DSMC simulations" => "overview_fixedweight.md",
            "Variable-weight DSMC simulations" => "overview_varweight.md",
            "1D DSMC simulations" => "overview_1d.md",
            "Fokker-Planck simulations" => "overview_fp.md",
        ],
        "Tutorials" => [
            "Contiguous indexing" => "contiguous_indexing.md",
            "Modelling ionization reactions" => "modelling_ionization.md",
            "Multithreaded simulations" => "multithreaded.md",
        ],
        "API reference" => [
            "Public API reference" => "reference_public.md",
            "Internal API reference" => "reference_internal.md",
        ],
        "Benchmarks" => "benchmarks.md"
    ]
)

deploydocs(
    repo = "github.com/merzbild/Merzbild.jl.git",
)