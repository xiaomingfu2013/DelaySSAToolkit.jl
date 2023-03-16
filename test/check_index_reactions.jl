using Catalyst
using Test
using DelaySSAToolkit
rn = @reaction_network begin
    @parameters ρ, r
    @species S, E, I, R
    ρ, S + I --> E + I
    1 / I, I --> R
    r, I --> R
end

m = get_reaction_idx(rn)
eqs = equations(rn)
@test m[:, 2] == [1, 3, 2]
@test m[:, 1] == eqs[[1, 3, 2]]
