module CompatUtils

using Core: UndefKeywordError
macro required(sym)
    sym
end

end  # module
