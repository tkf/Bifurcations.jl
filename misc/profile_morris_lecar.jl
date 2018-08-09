include("../benchmark/bench_morris_lecar.jl")

@time flc_solver = pre_solve!(make_morris_lecar_flc_solver())
@time sweep!(flc_solver)

@time flc_solver = pre_solve!(make_morris_lecar_flc_solver())
Profile.clear()
@profile sweep!(flc_solver)
