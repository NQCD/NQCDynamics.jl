import JSON
using CairoMakie
using Glob

# Get information from outside the plotting script: What is the newest version and what are we comparing to?
newest_version = get(ENV, "CURRENT_VERSION", "v0.15.1") # A known tested version
github_context = get(ENV, "NEW_VERSION", "new") # The new thing to test. 

# Base path of the benchmark repository
benchmark_path = get(ENV, "BENCHMARK_OUTPUTS_DIR")

