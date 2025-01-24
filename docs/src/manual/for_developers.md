# For developers

## Project structure

## Running tests locally

## Building the documentation locally

## Style guide

## Performance Benchmarks

There is a benchmark suite defined in `GalerkinToolkitExamples/benchmarks`. This uses [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) to perform the timings
and [github-action-benchmark](https://github.com/benchmark-action/github-action-benchmark) to collect the results and store them in the `benchmarks` branch. Upon merging of a PR, and after the docs have built, the benchmark results are copied from the `benchmarks` branch into the `gh-pages` branch. This is required to prevent overwriting of the previous results by docs deployment etc.

Graphs of performance changes over time (per commit hash) can then be viewed here: <https://galerkintoolkit.github.io/GalerkinToolkit.jl/dev/bench/>.

The github action can be configured (in `.github/workflows/CI.yml`, job `benchmark`) to fail if the performance change is beyond a given threshold. Look for the `alert-threshold:` and `fail-on-alert:` keys.

More benchmarks can be added (or existing ones modified) in `GalerkinToolkitExamples/benchmarks/run_benchmarks.jl`.
