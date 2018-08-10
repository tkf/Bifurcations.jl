module TestBifurcations
using Test
using Compat: @warn, stdout

RUNINFO = []

TEST_GROUPS = Dict{String, Vector{String}}(
    # ENV["CI_GROUP"] => test files
    "0" => [
        "test_utils.jl",
        "test_normal_form.jl",
        "test_smoke.jl",
        "test_jacobian.jl",
        "test_residual_jacobian.jl",
        "test_calcium.jl",
        "test_predator_prey.jl",
        "test_bazykin_85.jl",
        "test_bautin.jl",
        "test_examples.jl",
        "test_vs_svector.jl",
        "test_duffing_van_der_pol.jl",
        "test_hopf_to_lc.jl",
    ],
    "1" => [
        "test_reparametrization.jl",
        "test_reparametrized_bautin.jl",
    ],
    "2" => [
        "test_fold_lc.jl",
        "test_morris_lecar.jl",
    ],
)
TEST_GROUPS["all"] = vcat(last.(sort(collect(TEST_GROUPS), by=first))...)

@testset "$file" for file in TEST_GROUPS[get(ENV, "CI_GROUP", "all")]
    _, time, bytes, gctime, memallocs =  @timed @time include(file)
    push!(RUNINFO, Dict(
        :file => file,
        :time => time,
        :bytes => bytes,
        :gctime => gctime,
        :memallocs => memallocs,
    ))
end

function print_runinfo(io = stdout, runinfo = RUNINFO)
    io = IOContext(io, :compact => true)  # reduce number of digits shown
    titles = ["File", "Time (sec)", "%GC", "GiB"]
    table = Array{Any}((length(runinfo) + 1, length(titles)))
    table[1, :] .= titles
    for (i, rec) in enumerate(runinfo)
        table[i + 1, 1] = rec[:file]
        table[i + 1, 2] = rec[:time]
        table[i + 1, 3] = rec[:gctime] / rec[:time] * 100
        table[i + 1, 4] = rec[:bytes] / 2^30
    end
    Base.print_matrix(io, table)
    println(io)
end

try
    println()
    print_runinfo()
    println()
catch err
    @warn "Error while printing RUNINFO"
    @warn err
end

end  # module
