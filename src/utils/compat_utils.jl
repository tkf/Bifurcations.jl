module CompatUtils

@static if VERSION > v"0.7-"
    using Core: UndefKeywordError
    macro required(sym)
        sym
    end
else
    # UndefKeywordError added in: RFC: required keyword arguments
    # https://github.com/JuliaLang/julia/pull/25830/files
    struct UndefKeywordError <: Exception
        var::Symbol
    end
    Base.showerror(io::IO, ex::UndefKeywordError) =
        print(io, "UndefKeywordError: keyword argument $(ex.var) not assigned")

    macro required(sym)
        Expr(:kw, esc(sym), :(throw(UndefKeywordError($(QuoteNode(sym))))))
    end
end

end  # module
