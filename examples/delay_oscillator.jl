using DelaySSAToolkit
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

delay_trigger_affect! = function (de_chan, rng)
   append!(de_chan[1], 10.)
end

delay_trigger =Dict(2=>delay_trigger_affect!) # 可以 trigger 多个delay_reactions 这个表示 10s 后同时增加作用 两次 1 号delay 反应 
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
tf = 20.
de_chan0 = [[10.]]
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
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpsets, delaysets, de_chan0, save_positions = (false,false))
# jprob.discrete_jump_aggregation
sol = solve(djprob, SSAStepper(),seed=1234, saveat =.02)

djprob2 = DelayJumpProblem(dprob, DelayRejection(), jumpsets, delaysets, de_chan0, save_positions = (true,true))
sol2 = solve(djprob2, SSAStepper(),seed=1234)


plot(sol.t,[sol.u[i][1] for i in 1:length(sol.u)])
scatter!(sol2.t,[sol2.u[i][1] for i in 1:length(sol2.u)],legend = false)
plot!(sol.t,[sol.u[i][2] for i in 1:length(sol.u)])
scatter!(sol2.t,[sol2.u[i][2] for i in 1:length(sol2.u)],legend = false)
