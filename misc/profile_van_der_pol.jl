include("../benchmark/bench_van_der_pol.jl")

@time lc_solver = pre_solve!(make_van_der_pol_lc_solver())
@time sweep!(lc_solver)

@time lc_solver = pre_solve!(make_van_der_pol_lc_solver())
Profile.clear()
@profile sweep!(lc_solver)
