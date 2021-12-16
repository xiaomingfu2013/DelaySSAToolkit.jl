"""
-`delay_trigger::Dict{Int,Any}`     #  those reactions in the Markovian part that trigger the change of the state of the delay channels or/and the state of the reactants upon initiation. Keys: Indices of reactions defined in the Markovian part that can trigger the delay reactions; Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

-`delay_interrupt::Dict{Int,Any}`   # those reactions in the Markovian part that change the state of the delay channels or/and the state of the reactants in the middle of on-going delay reactions. Keys: Indices of reactions defined in the Markovian part that can interrupt the delay reactions; Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

-`delay_complete::Dict{Int,Any}`    # those reactions that are initiated by delay trigger reactions and change the state of the delay channels or/and the state of the reactants upon completion. Keys: Indices of the delay channel; Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

-`delay_trigger_set::Vector{Int}`   # collect(keys(delay_trigger))

-`delay_interrupt_set::Vector{Int}` # collect(keys(delay_interrupt))
"""
mutable struct DelayJumpSet
    delay_trigger::Dict{Int,Any}
    delay_complete::Dict{Int,Any}
    delay_interrupt::Dict{Int,Any}
    delay_trigger_set::Vector{Int}
    delay_interrupt_set::Vector{Int}
end
DelayJumpSet(delay_trigger,delay_complete,delay_interrupt) = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt, collect(keys(delay_trigger)), collect(keys(delay_interrupt)))

