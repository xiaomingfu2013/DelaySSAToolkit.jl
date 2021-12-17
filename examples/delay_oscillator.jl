using DelaySSAToolkit
using Catalyst
# using Plots
################################################################################
# the following example illustrates how to deal with a oscillatory system
# d: X-> 0; k1/(1+Y^2): 0-> X; [trigger X-> Y after τ time;] k2*Y/(1+Y):  Y-> 0;  

rn = @reaction_network begin
    1/(1+Y^2), 0 --> X
    1/(1+Y),   Y --> 0
end
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)

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
delay_complete = Dict(1=>[2=>1, 1=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(true,true))
sol = solve(djprob, SSAStepper(),seed=1234)

djprob2 = DelayJumpProblem(jumpsys, dprob, DelayDirect(), delayjumpset, de_chan0, save_positions = (true,true))
sol2 = solve(djprob2, SSAStepper(),seed=1234)

plot(sol.t,[sol.u[i][1] for i in 1:length(sol.u)])
plot!(sol.t,[sol.u[i][2] for i in 1:length(sol.u)])
scatter(sol2.t,[sol2.u[i][1] for i in 1:length(sol2.u)],legend = false)
scatter(sol2.t,[sol2.u[i][2] for i in 1:length(sol2.u)],legend = false)
