using BenchmarkTools
const SUITE = BenchmarkGroup()
SUITE["morris_lecar"] = include("bench_morris_lecar.jl")
