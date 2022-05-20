using DelaySSAToolkit, Catalyst

rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r

jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
states(jumpsys)

u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)


ps = [1e-4, 1e-2]
τ = 20.
dprob = DiscreteProblem(jumpsys,u0,tspan, ps)
delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
# delay_trigger = Dict(1=>[1=>τ])
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
# algo = DelayDirect()
# algo = DelayDirectCR()
algo = DelayRejection()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
jprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0, save_positions=(true,true))
@time sol = solve(jprob, SSAStepper(), seed = 1234)
using Plots; theme(:vibrant)
fig = plot(sol, label = ["S" "I" "E" "R"], linewidth = 3, legend = :top, ylabel = "# of individuals", xlabel = "Time", fmt=:png)
savefig(fig,"docs/src/assets/seir.png")