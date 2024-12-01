# Developers guide

## Performance Benchmarks
There is a benchmark suite defined in `GalerkinToolkitExamples/benchmarks`. This uses [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) to perform the timings
and [github-action-benchmark](https://github.com/benchmark-action/github-action-benchmark) to collect the results and store them in the `gh-pages` branch. Graphs of performance
changes over time (per commit hash) can be viewed here: <https://galerkintoolkit.github.io/GalerkinToolkit.jl/dev/bench/>.

More benchmarks can be added (or existing ones modified) in `GalerkinToolkitExamples/benchmarks/run_benchmarks.jl`.
