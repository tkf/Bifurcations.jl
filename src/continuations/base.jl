abstract type AbstractContinuationProblem{iip} end
abstract type AbstractProblemCache{P} end

function get_prob_cache end
function get_u0 end
function residual end
function residual! end
function residual_jacobian! end
function isindomain end
