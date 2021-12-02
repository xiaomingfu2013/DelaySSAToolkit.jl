using DelaySSAToolkit, DiffEqJump
# Gon, Goff, N, M, P
#     k10, Gon + P --> Goff
#     k01, Goff  --> Gon + P
#     ρ1, Gon  --> Gon ==>τ1 N, N ==>τ2 M
#     ρ2, M  --> M + P
#     d1, N --> ∅
#     d2, M --> ∅
#     d3, P --> ∅  
#     
rate = [1.,0.1,5.,1.,0.01,0.1,0.1]
reactant_stoich = [[1=>1,5=>1],[2=>1],[1=>1],[4=>1],[3=>1],[4=>1],[5=>1]]
net_stoich = [[1=>-1,2=>1,5=>-1],[1=>1,2=>-1,5=>1],[],[5=>1],[3=>-1],[4=>-1],[5=>-1]]
mass_jump = MassActionJump(rate, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,mass_jump)

# reaction_idx =>  Dict( delay_react_idx => delay_time )
τ1, τ2 = [1., 10.]
delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], τ1)
end
# delay_trigger =Dict(3=>delay_trigger_affect!) # 可以 trigger 多个delay_reactions 这个表示 10s 后同时增加作用 两次 1 号delay 反应
delay_trigger =Dict(3=>[1=>τ1]) # 可以 trigger 多个delay_reactions 这个表示 10s 后同时增加作用 两次 1 号delay 反应
delay_complete_affect1! = function (integrator, rng)
    integrator.u[3] += 1
    append!(integrator.de_chan[2], τ2)
end
delay_complete = Dict(1=>delay_complete_affect1!, 2=>[3=>-1,4=>1])
delay_interrupt_affect! = function (integrator, rng)
   i = rand(rng, 1:length(integrator.de_chan[2])) # second channel is for N => M
   deleteat!(integrator.de_chan[2],i)
end
delay_interrupt = Dict(5=>delay_interrupt_affect!) # reactions 能够对 delay channel中造成影响 keys: reaction idx -> values: delay_affect! 
delaysets = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)
using Random
u0 = [1,0,2,0,0]
tf = 500.
de_chan0 = [[0.1,0.2,0.3],[5.,9.]]  #first being production of N delay channel, second being N ==> M delay channel 
# de_chan0 = [[],[]]  #first being production of N delay channel, second being N ==> M delay channel 
tspan = (0.,tf)
dprob = DiscreteProblem(u0, tspan)
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delaysets, de_chan0, save_positions = (false,false))
delay_sol =@time solve(djprob, SSAStepper(), seed=2, saveat=1., save_delay_channel=true)
delay_sol.chan_sol

show(delay_sol)

delay_sol.chan_sol


n_N = delay_sol[3,:]
n_M = delay_sol[4,:]
n_P = delay_sol[5,:]
ens_prob = EnsembleProblem(djprob)
@time ens = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=10^4)


supertype(typeof(delay_sol))

using Plots; theme(:vibrant)
plot(0:1:500,[n_N n_M n_P],label=["N" "M" "P"], xlabel = "Time", ylabel = "# of reactants")
using DifferentialEquations.EnsembleAnalysis
slice_end = componentwise_vectors_timepoint(ens, tf)
histogram(slice_end[5], bins = 0:1:100)

