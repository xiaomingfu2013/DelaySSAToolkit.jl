macro inherit(name, base, fields)
    base_type = Core.eval(@__MODULE__, base)
    base_fieldnames = fieldnames(base_type)
    base_types = [t for t in base_type.types]
    base_fields = [:($f::$T) for (f, T) in zip(base_fieldnames, base_types)]
    res = :(mutable struct $name end)
    push!(res.args[end].args, base_fields...)
    push!(res.args[end].args, fields.args...)
    return res
end
using DiffEqBase, DiffEqJump
u0 = [1]
tspan = (0,1.)
p = 0.
prob = DiscreteProblem(u0, tspan, p)
alg = SSAStepper()
t = [prob.tspan[1]]
u = [copy(prob.u0)]

sol = DiffEqBase.build_solution(prob,alg,t,u,dense=false,
                         calculate_error = false,
                         destats = DiffEqBase.DEStats(0),
                         interp = DiffEqBase.ConstantInterpolation(t,u))
sol|>typeof

@inherit DSSASolution{deType} ODESolution begin
    de_chan::deType
end

# ODESolution|>typeof

Core.eval(@__MODULE__, ODESolution)
ODESolution


delay_trigger_aff = [1=>1.,2=>3.]

for (chan_idx, τ) in delay_trigger_aff
    println([chan_idx], [τ])
end