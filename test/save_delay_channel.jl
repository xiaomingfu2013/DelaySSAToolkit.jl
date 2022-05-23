using DelaySSAToolkit
using Test, DiffEqJump, StaticArrays
#σ_off: Gon -> Goff
#σ_on: Goff -> Gon
#ρ_on: Gon -> Gon + N, triggers N => τ 0
#d: N -> 0
# 1. Gon, 2. Goff, 3. N

params = [0.1, 0.1, 10.0, 0.1, 10.0, 50.0]
σ_off, σ_on, ρ_on, d, τ, tf = params
rates = [σ_off, σ_on, ρ_on, d]
react_stoich = [[1 => 1], [2 => 1], [1 => 1], [3 => 1]]
net_stoich = [[1 => -1, 2 => 1], [1 => 1, 2 => -1], [3 => 1], [3 => -1]]
mass_action_jump = MassActionJump(rates, react_stoich, net_stoich; scale_rates=false)
jumpset = JumpSet((), (), nothing, mass_action_jump)
delay_trigger = Dict(3 => [1 => τ])
delay_complete = Dict(1 => [3 => -1])
delay_affect! = function (integrator, rng)
    i = rand(rng, 1:length(integrator.de_chan[1]))
    deleteat!(integrator.de_chan[1], i)
end
delay_interrupt = Dict(4 => delay_affect!)

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
u0 = [0, 1, 0]
de_chan0 = [[]]
tspan = (0.0, tf)
dprob = DiscreteProblem(u0, tspan)
algs = [DelayRejection(), DelayDirect(), DelayMNRM(), DelayDirectCR()]


## TEST
# alg = algs[3]
# djprob = DelayJumpProblem(dprob, alg, jumpset, delayjumpset, de_chan0, save_positions=(false, false), save_delay_channel=true)

# dep_gr = djprob.jump_callback.discrete_callbacks[1].affect!.dep_gr
# p = djprob.jump_callback.discrete_callbacks[1].affect!
# DelaySSAToolkit.var_to_jumps_map(3, mass_action_jump)
# integrator = DelaySSAToolkit.init(djprob, SSAStepper())
# DelaySSAToolkit.vars_to_jumps_delay(p, integrator, [3])
# DelaySSAToolkit.dep_gr_delay(p, integrator)
# sol = solve(djprob, SSAStepper(), seed = 2)
# for i in eachindex(sol.u)
#     @test sol.u[i][3] == length(sol.channel[i][1])
# end  

for alg in algs
    djprob = DelayJumpProblem(dprob, alg, jumpset, delayjumpset, de_chan0, save_positions=(true, true), save_delay_channel=true)
    @info "Testing method $(alg)"
    sol = solve(djprob, SSAStepper())
    
    for i in eachindex(sol.u)
        @test sol.u[i][3] == length(sol.channel[i][1])
    end
end

