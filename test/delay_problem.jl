using Test, DelaySSAToolkit, Catalyst

ρ, r = [1e-4, 1e-2]
rates = [ρ, r]
reactant_stoich = [[1=>1,2=>1],[2=>1]]
net_stoich = [[1=>-1,3=>1],[2=>-1,4=>1]]
mass_jump = MassActionJump(rates, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,mass_jump)

u0 = [999,1,0,0] # S, I, E, R
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2] # parameters for ρ, r
τ = 20.
dprob = DiscreteProblem(u0,tspan,ps)

delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
de_chan0 = [[]]
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions=(false, false), save_delay_channel = true)
sol = solve(djprob, SSAStepper(), seed = 1, saveat =10.)

rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
djprob_ = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(false,false), save_delay_channel = true)

sol_ = solve(djprob_, SSAStepper(), seed =1, saveat = 10.)

@testset for i in 1:length(sol.t)
    @test sol[i] == sol_[i]
end

@testset for i in 1:3
    @test sol[i,:] == sol_[i,:]
    @test sol(:u, i) == sol_(:u, i)
end
@testset for t in sol.t
    @test sol(t) == sol_(t)
end
@test sol(:channel, 1) == sol_(:channel, 1)