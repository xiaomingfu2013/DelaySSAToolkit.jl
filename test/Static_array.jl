using Catalyst
using StaticArrays
using DelaySSAToolkit
# using DifferentialEquations

rn = @reaction_network begin
    1, 0 --> A
end

delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], 10)
end

u0 =@SVector [ 0 ]
tspan = (0., 200.)

delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[1=>-1])
delay_interrupt = Dict() 
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

dprob = DiscreteProblem(rn, u0, tspan)
jsys = convert(JumpSystem, rn)
djprob = DelayJumpProblem(jsys, dprob, DelayRejection(), delayjumpset, [[]])
sol = solve(djprob, SSAStepper())
using Plots
plot(sol)