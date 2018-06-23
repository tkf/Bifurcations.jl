tkindstr(x) = tkindstr(timekind(x))
tkindstr(::Discrete) = "Discrete"
tkindstr(::Continuous) = "Continuous"

function print_header(io::IO, point::Union{SpecialPoint,
                                           SpecialPointInterval})
    print(io, nameof(typeof(point)), " <",
          tkindstr(point), " ",
          point.point_type,
          ">")
end

set_if_not(io, key, val) = haskey(io, key) ? io : IOContext(io, key => val)

function Base.show(io::IO, point::SpecialPoint)
    print_header(io, point)
    println(io)
    if ! get(io, :compact, false)
        io = set_if_not(io, :compact, true)  # reduce number of digits shown
        println(io, "u = ", point.u)
    end
end

function Base.show(io::IO, point::SpecialPointInterval)
    print_header(io, point)
    println(io)
    if ! get(io, :compact, false)
        io = set_if_not(io, :compact, true)  # reduce number of digits shown
        println(io, "happened between:")
        println(io, "  u0 = ", point.u0)
        println(io, "  u1 = ", point.u1)
    end
end
