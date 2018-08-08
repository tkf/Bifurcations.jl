using BenchmarkTools
const SUITE = BenchmarkGroup()
include("bench_van_der_pol.jl")
include("bench_morris_lecar.jl")
