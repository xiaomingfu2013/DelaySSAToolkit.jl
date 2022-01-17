using DelaySSAToolkit, Catalyst
rn = @reaction_network begin
    1/(1+Y^2), 0 --> X
    1/(1+Y),   Y --> 0
end
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)
# states(rn)
u0 = [0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
τ = 20.
dprob = DiscreteProblem(jumpsys, u0, tspan)
# jumpsys.dep_graph

delay_trigger_affect! = function (integrator, rng)
  append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 1=>-1], 2=>[2=>1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
# alg = DelayDirectCR()
alg = DelayMNRM()
# alg = DelayRejection()
# alg = DelayDirect()
djprob = DelayJumpProblem(jumpsys, dprob, alg, delayjumpset, de_chan0)

p = djprob.discrete_jump_aggregation;
# @test p.dep_gr = 
integrator = DelaySSAToolkit.init(djprob, SSAStepper())
integrator.de_chan
DelaySSAToolkit.dep_gr_delay(p, integrator)

