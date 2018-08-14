module TestMorrisLecar
include("preamble.jl")

using StaticArrays: SVector

using Bifurcations: resolved_points
using Bifurcations.BifurcationsBase: count_special_points
using Bifurcations.Codim2LimitCycle: FoldLimitCycleProblem
using Bifurcations.Examples.MorrisLecar

# To test start_from_nearest_root, let's start from an arbitrary
# place:
prob = MorrisLecar.make_prob(
    MorrisLecar.MorrisLecarParam(y=0.1, z=0.1);
    u0 = SVector(0.0, 0.0),
)
solver1 = init(
    prob;
    start_from_nearest_root = true,
)
@time solve!(solver1)

all_specials = special_points(solver1)
all_hopf = special_points(solver1, Codim1.PointTypes.hopf)
all_sn = special_points(solver1, Codim1.PointTypes.saddle_node)
@test length(all_specials) == 3
@test length(all_sn) == 2
@test length(all_hopf) == 1

hopf_point, = all_hopf
hopf_point

codim2_prob = BifurcationProblem(
    hopf_point,
    solver1,
    (@lens _.z),
    (-1.0, 1.0),
)
hopf_solver = init(
    codim2_prob;
    nominal_angle_rad = 0.01,
)
solve!(hopf_solver)

all_bautin = special_points(hopf_solver, Codim2.PointTypes.bautin)
@assert length(all_bautin) == 1
bautin_point, = all_bautin
bautin_point

using Bifurcations.Codim2LimitCycle: FoldLimitCycleProblem

@time flc_prob = FoldLimitCycleProblem(
    bautin_point,
    hopf_solver;
    num_mesh = 30,
    degree = 4,
);
flc_solver = init(
    flc_prob;
    start_from_nearest_root = true,
    max_branches = 0,
    max_samples = 20,
    bidirectional_first_sweep = false,
    nominal_angle_rad = 2π * (5 / 360),
    verbose = true,
)
@time solve!(flc_solver)



sn_prob = BifurcationProblem(
    special_points(solver1, Codim1.PointTypes.saddle_node)[1],
    solver1,
    (@lens _.z),
    (-1.0, 1.0),
)
sn_solver = init(
    sn_prob;
    nominal_angle_rad = 2π * (5 / 360),
    max_samples = 1000,
    start_from_nearest_root = true,
)
@time solve!(sn_solver)

@test count_special_points(sn_solver) == Dict(
    Codim2.PointTypes.bogdanov_takens => 1,
    Codim2.PointTypes.cusp => 1,
)


hopf_prob2 = BifurcationProblem(
    special_points(sn_solver, Codim2.PointTypes.bogdanov_takens)[1],
    sn_solver,
)
hopf_solver2 = init(hopf_prob2)
@time solve!(hopf_solver2)

count_hopf_solver2 = count_special_points(hopf_solver2)

@test Set(keys(count_hopf_solver2)) <= Set([Codim2.PointTypes.bogdanov_takens])
if Codim2.PointTypes.bogdanov_takens in keys(count_hopf_solver2)
    @test count_hopf_solver2 == Dict(
        Codim2.PointTypes.bogdanov_takens => 1,
    )
end

end  # module
