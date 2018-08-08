using BenchmarkTools
const SUITE = BenchmarkGroup()
SUITE["van_der_pol"] = include("bench_van_der_pol.jl")
SUITE["morris_lecar"] = include("bench_morris_lecar.jl")
