using DelaySSAToolkit, Catalyst
using Test
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


delay_trigger_affect! = function (integrator, rng)
  append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1], 2=>[1=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
# alg = DelayDirectCR()
# alg = DelayMNRM()
for alg = [DelayDirectCR(), DelayMNRM()]
  djprob = DelayJumpProblem(jumpsys, dprob, alg, delayjumpset, de_chan0)

  p = djprob.discrete_jump_aggregation;
  integrator = DelaySSAToolkit.init(djprob, SSAStepper())
  # integrator.de_chan
  # delay_complete[1]
  # p.dep_gr
  # DelaySSAToolkit.dep_gr_delay(p, integrator)
  
  @test p.vartojumps_map == [[],[1,2]]
  @test p.dep_gr == [[1],[1,2]]
  @test DelaySSAToolkit.dep_gr_delay(delayjumpset, p.vartojumps_map, 2) == Dict(1=>[1,2], 2=>[])
end