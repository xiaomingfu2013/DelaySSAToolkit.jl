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
u0 =@SVector [0, 1, 0]
de_chan0 = [[]]
tspan = (0.0, tf)
dprob = DiscreteProblem(u0, tspan)
algs = [DelayRejection(), DelayDirect(), DelayMNRM(), DelayDirectCR()]
alg = algs[2]
djprob = DelayJumpProblem(dprob, alg, jumpset, delayjumpset, de_chan0, save_positions=(false, false), save_delay_channel=true)
# dep_gr = djprob.jump_callback.discrete_callbacks[1].affect!.dep_gr

sol = solve(djprob, SSAStepper(), seed = 1)

sol.u[end]
sol.channel[end][1]

for i in eachindex(sol.u)
    @test sol.u[i][3] == length(sol.channel[i][1])
end  

for alg in algs
    djprob = DelayJumpProblem(dprob, alg, jumpset, delayjumpset, de_chan0, save_positions=(true, true), save_delay_channel=true)
    @info "Testing method $(alg)"
    sol = solve(djprob, SSAStepper(), seed = 1)
    
    for i in eachindex(sol.u)
        @test sol.u[i][3] == length(sol.channel[i][1])
    end
end


u0_ = [0,1,0]
dprob_ = DiscreteProblem(u0_, tspan)
rates_ = [σ_off, σ_on, ρ_on]
react_stoich_ = [[1 => 1], [2 => 1], [1 => 1]]
net_stoich_ = [[1 => -1, 2 => 1], [1 => 1, 2 => -1], [3 => 1]]
mass_action_jump_ = MassActionJump(rates_, react_stoich_, net_stoich_; scale_rates=false)
affect! = function (integrator)
    integrator.u[3] -= 1
end
rate2 =  (u,p,t) -> 0.1*u[3]
constant_rate_jump = ConstantRateJump(rate2, affect!)
jumpset_ = JumpSet((), (constant_rate_jump,), nothing, mass_action_jump_)



djprob_ = DelayJumpProblem(dprob_, DelayRejection(), jumpset_, delayjumpset, de_chan0, save_positions=(true, true), save_delay_channel=true)

sol = solve(djprob_, SSAStepper(), seed= 1)
sol.u
sol.channel

dep_gr = [[1,2,3],[1,2,3],[3,4],[4]]
for alg in algs
    djprob = DelayJumpProblem(dprob_, alg, jumpset_, delayjumpset, de_chan0, save_positions=(true, true), save_delay_channel=true, dep_graph = dep_gr)
    @info "Testing method $(alg)"
    sol = solve(djprob, SSAStepper(), seed = 1)
    
    for i in eachindex(sol.u)
        @test sol.u[i][3] == length(sol.channel[i][1])
    end
end