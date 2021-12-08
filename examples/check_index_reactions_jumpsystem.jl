using Catalyst
rn = @reaction_network begin
    ρ, S+I --> E+I
    1/I, I --> R
    r, I --> R
    # t, I --> R
end ρ r
rn.eqs
jumpsys.eqs

jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
rn.eqs

jumpsys.eqs

equations(jumpsys).x[1]
equations(jumpsys).x[2]

equations(jumpsys).x[3]


function get_reaction_idx(rn::Catalyst.ReactionSystem)
    massaction_order = []
    constantrate_order = []
    rxvars = []
    for (i,rx) in enumerate(ModelingToolkit.get_eqs(rn))
        (rx.rate isa ModelingToolkit.Symbolic) && ModelingToolkit.get_variables!(rxvars, rx.rate)
        @inbounds for i = 1:length(rxvars)
            if isequal(rxvars[i], ModelingToolkit.get_iv(rn))
                error("Does not support VariableRateJump")
                break
            end
        end
        if Catalyst.ismassaction(rx, rn)
            push!(massaction_order, i)
        else
            push!(constantrate_order, i)
        end
    end
    order = vcat(massaction_order,constantrate_order)
    vcat( [[rn.eqs[i] order[i]] for i in eachindex(order)]...) 
end
get_reaction_idx(rn)
