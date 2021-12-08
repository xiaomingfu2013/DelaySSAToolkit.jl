using Catalyst
rn = @reaction_network begin
    ρ, S+I --> E+I
    1/I, I --> R
    r, I --> R
    # t, I --> R
end ρ r


jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
rn.eqs
jumpsys.eqs

equations(jumpsys).x[1]
equations(jumpsys).x[2]

equations(jumpsys).x[3]

"""
    function get_reaction_idx(rn::Catalyst.ReactionSystem)

Get the rearranged order for a reaction system, returns a matrix `m` with first column as the rearranged reactions, second column is the order in the orign `ReactionSystem`
If you know the old idx for a certain reaction in `ReactionSystem`, you can use 
`m[:,2][idx]` to get the new idx. 
"""
function get_reaction_idx(rn::Catalyst.ReactionSystem)
    massactionjump_order = Int64[]
    constantjump_order = Int64[]
    rxvars = []
    for (i,rx) in enumerate(ModelingToolkit.get_eqs(rn))
        empty!(rxvars)
        (rx.rate isa ModelingToolkit.Symbolic) && ModelingToolkit.get_variables!(rxvars, rx.rate)
        @inbounds for j = 1:length(rxvars)
            if isequal(rxvars[j], ModelingToolkit.get_iv(rn))
                error("Does not support VariableRateJump")
            end
        end
        if Catalyst.ismassaction(rx, rn)
            push!(massactionjump_order, i)
        else
            push!(constantjump_order, i)
        end
    end
    order = vcat(massactionjump_order, constantjump_order)
    # println("The rearranged order :")
    vcat([[rn.eqs[order[i]] order[i]] for i in eachindex(order)]...) 
end


m = get_reaction_idx(rn)
m[:,2][3]

