# Merzbild.jl AGENTS.md file

## General instructions and tips

* Ignore all directories and files listed in `.gitignore`.
* Project documentation apart from docstrings is in `docs/src`.
* Particle indexing documentation is split across `docs/src/overview_blocks` and `docs/src/contiguous_indexing.md`.
* Memory allocations are to be avoided at all costs (unless new instances of structs are instantiated or completely new particles add),
prefer explicit loops and re-using data to cleaner-looking operations using vectorized syntax.
* NEVER modify `Project.toml` and `Manifest.toml` by yourself
* Only include comments when absolutely necessary. When the function name or implementation clearly indicates its purpose or behavior, redundant comments are unnecessary.

## Testing instructions

* Run `julia --project=. -e 'using Pkg; Pkg.test()'`to run all project tests.
