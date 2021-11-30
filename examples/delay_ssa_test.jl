using DelaySSA
# using Plots
################################################################################
# the following example illustrates how to deal with a oscillatory system
# d: X-> 0; k1/(1+Y^2): 0-> X; [trigger X->τ 2*Y;] k2*Y/(1+Y):  Y-> 0;  


rate1 = [1e-2]
reactant_stoch = [[1=>1]]
net_stoch = [[1=>-1]]
mass_jump = MassActionJump(rate1, reactant_stoch, net_stoch; scale_rates =false)

rate2 = (u,p,t) -> 1/(1+u[2]^2)
affect2! = function (integrator)
  integrator.u[1] += 1
end
cons_jump1 = ConstantRateJump(rate2,affect2!)

rate3 = (u,p,t) -> u[2]/(1+u[2])
affect3! = function (integrator)
  integrator.u[2] -= 1
end
cons_jump2 = ConstantRateJump(rate3,affect3!)

jumpsets = JumpSet((),(cons_jump1,cons_jump2,),nothing,[mass_jump])

# DelaySSA.get_num_majumps(jumpsets)
# reaction_idx =>  Dict( delay_react_idx => delay_time )
delay_trigger =Dict(2=>Dict(1=>[10.])) # 可以 trigger 多个delay_reactions 这个表示 10s 后同时增加作用 两次 1 号delay 反应 
delay_complete = Dict(1=>[1=>-1,2=>1])

# for (i,ξ) in delay_complete[1]
#     println(i,ξ)
# end
# length(Dict(2=>1,3=>4))

delay_affect! = function (de_chan, rng)
   i = rand(rng, 1:length(de_chan[1]))
   deleteat!(de_chan[1],i)
end
# typeof(delay_affect!)
delay_period = Dict(1=>delay_affect!) # reactions 能够对 delay channel中造成影响的 keys : reaction idx -> delay idx + how many molecules will be consumed 
# delay_period[1](Dict(1=>[1,2,3]))

delaysets = DelayJumpSet(delay_trigger,delay_complete,delay_period)

u0 = [1,0]
tf = 15.
saveat = 1.
timepoints = Array(0:saveat:tf)
de_chan0 = Dict(1 =>[10.])
delaysys = DelayJumpSystem(jumpsets,delaysets)
delaysys.delayjumpsets


p = 0.
tspan = (0.,tf)
aggregatoralgo = DelayRejection()
save_positions = (false,false)
# using Random
# integrator = DelaySSA.initialize_delay_integrator(u0,de_chan0,delaysys,p,tspan,saveat)
# aggregator = DelaySSA.initialize_delay_aggregator(integrator,aggregatoralgo,delaysys,save_positions,Random.seed!(1234))

# aggregator(integrator)
# integrator
# aggregator

dprob = DiscreteProblem(u0, tspan, p)
jprob = JumpProblem(dprob, DelayRejection(), jumpsets, save_positions = (false,false))
# jprob.discrete_jump_aggregation

# jprob.u0
# include("../src/delaySSA_stepper.jl")

# de_chan0
djprob = DelayJumpProblem(jprob,delaysets,de_chan0)
# djprob.jump_callback.discrete_callbacks[end]
# # DelaySSA.DSSASolution([tspan[1]],[u0],[de_chan0])
# # jumpsets.variable_jumps
# djprob.delayjumpsets
# integrator = DiffEqBase.__init(djprob,SSAStepper(),seed = 1234)

# aggregator = integrator.cb.affect!
# aggregator.next_jump_time
# aggregator(integrator)
# integrator
# aggregator
# integrator.tstops

# # DiffEqBase.step!(integrator)
# integrator
# typeof(integrator.opts.callback.discrete_callbacks)<:Tuple{}


# DelaySSA.add_tstop!(integrator,10.)

# de_chan0
# DelaySSA.initialize_delay_integrator(u0,de_chan0,delaysys,p,tspan,1.)

sol, chan_sol = solve(djprob, SSAStepper(),seed=1234, saveat =1.)
sol.u[11]
chan_sol[11][1]
# sol[11]
# chan_sol[11][1]
# sol(30)
sol.u
chan_sol

delay_sol =@time DelaySSAsolve(delaysys,aggregatoralgo,u0,de_chan0,p,(0.,10.),seed=1234,save_delay_chan=true,saveat=.1)
delay_sol_nosaveat =@time DelaySSAsolve(delaysys,aggregatoralgo,u0,de_chan0,p,(0.,10.),seed=1234,save_delay_chan=true,saveat=nothing)
# 解决如何 更新 save_at
# de_chan0


delay_sol.t
delay_sol.u
delay_sol.de_chan

delay_sol_nosaveat.t
delay_sol_nosaveat.u
delay_sol_nosaveat.de_chan

plot(delay_sol.t,[delay_sol.u[i][1] for i in 1:length(delay_sol.u)])
scatter!(delay_sol_nosaveat.t,[delay_sol_nosaveat.u[i][1] for i in 1:length(delay_sol_nosaveat.u)],legend = false)
# copy(de_chan0)

sim =@time DelaySSA_multithreads(delaysys,aggregatoralgo,u0,de_chan0,p,(0.,100.),1e4,saveat=1.)
sim[1].u
sol1 = [[sim[i].u[t][1] for i in Int.(1:1e4)] for t in 1:1:101]
sol2 = [[sim[i].u[t][2] for i in Int.(1:1e4)] for t in 1:1:101]

plot(0:1:100,mean.(sol1))
plot!(0:1:100,mean.(sol2))