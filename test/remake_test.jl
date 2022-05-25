using Test, DelaySSAToolkit, Catalyst

rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r
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
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(false,false), save_delay_channel = true)

ps_ = 2*ps
de_chan0_ = [[1]]
u0_ = [998,1,1,0]
tspan_ = (0.,200.)

delay_trigger_ = Dict(1=>[1=>20.])
delay_complete_ = Dict(1=>[2=>2, 3=>-2])
# delay_interrupt_ = Dict()
djprob_ = remake(djprob, p = ps_, de_chan0 = de_chan0_, u0= u0_, tspan = tspan_, delay_trigger = delay_trigger_, delay_complete = delay_complete_)

@test djprob_.prob.p == ps_
@test djprob_.prob.u0 == u0_
@test djprob_.prob.tspan == tspan_
@test djprob_.delayjumpsets.delay_trigger == delay_trigger_
@test djprob_.delayjumpsets.delay_complete == delay_complete_
@test djprob_.de_chan0 == de_chan0_