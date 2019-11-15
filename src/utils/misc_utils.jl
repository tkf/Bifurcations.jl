module MiscUtils

import ProgressLogging

"""
    @progress_if cond [options...] expr

If `cond` evaluates to `true`, run `@progress(options..., expr)`; otherwise
just run `expr`.

`options` and `expr` are passed to
[`ProgressLogging.@progress`](https://junolab.org/ProgressLogging.jl/dev/#ProgressLogging.@progress-Tuple)
"""
macro progress_if(cond, args...)
    quote
        if $cond
            $ProgressLogging.@progress($(args...))
        else
            $(args[end])
        end
    end |> esc
end

end  # module
