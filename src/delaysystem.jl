"""
$(TYPEDEF)

A delay jump set...

# Fields
$(FIELDS)

# Notes
- `delay_trigger::Dict{Int,Any}`:    those reactions in the Markovian part that trigger the change of the state of the delay channels or/and the state of the reactants upon initiation. 
    - Keys: Indices of reactions defined in the Markovian part that can trigger the delay reactions; 
    - Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

- `delay_interrupt::Dict{Int,Any}`:  those reactions in the Markovian part that change the state of the delay channels or/and the state of the reactants in the middle of on-going delay reactions. 
    - Keys: Indices of reactions defined in the Markovian part that can interrupt the delay reactions; 
    - Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

- `delay_complete::Dict{Int,Any}`:     those reactions that are initiated by delay trigger reactions and change the state of the delay channels or/and the state of the reactants upon completion. 
    - Keys: Indices of the delay channel; 
    - Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

- `delay_trigger_set::Vector{Int}`:  collect(keys(delay_trigger))

- `delay_interrupt_set::Vector{Int}`:  collect(keys(delay_interrupt))
"""
mutable struct DelayJumpSet
    """those reactions in the Markovian part that trigger the change of the state of the delay channels or/and the state of the reactants upon initiation."""
    delay_trigger::Dict{Int,Any}
    """those reactions in the Markovian part that change the state of the delay channels or/and the state of the reactants in the middle of on-going delay reactions."""
    delay_complete::Dict{Int,Any}
    """those reactions that are initiated by delay trigger reactions and change the state of the delay channels or/and the state of the reactants upon completion."""
    delay_interrupt::Dict{Int,Any}
    """collect(keys(delay_trigger))"""
    delay_trigger_set::Vector{Int}
    """collect(keys(delay_interrupt))"""
    delay_interrupt_set::Vector{Int}    
end

DelayJumpSet(delay_trigger,delay_complete,delay_interrupt) = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt, collect(keys(delay_trigger)), collect(keys(delay_interrupt)))




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

export get_reaction_idx