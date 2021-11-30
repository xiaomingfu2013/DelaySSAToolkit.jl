using DelaySSAToolkit, Catalyst
using DiffEqJump

rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r

jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2]
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
τ = 20.
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
jprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(true,true))
sol = solve(jprob, SSAStepper())
fig = plot(sol, label = ["S" "I" "E" "R"], linewidth = 3, legend = :top, ylabel = "Cases", xlabel = "Time", fmt=:svg)
savefig(fig,"docs/src/assets/seir.svg")