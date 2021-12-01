# The model is 
# C: 0 -> X_A, 
# γ: X_A -> 0 
# β: X_A ->X_I1+X_I2 : trigger X_I1, X_I2 -> 0 after τ time and 
# γ: X_I1 ->0 : delay_interrupt
# γ: X_I2 ->0 : delay_interrupt


# According to Lafuerza, Toral 2011 the exact solution is obtainable
# if we denote x_A(t) to be the mean value of X_A at time t, and x_I(t) the mean value of X_I at time t, then 
# x_A(t)= C/a (1-e^(-at))  
# x_I(t)= 0<=t<=τ ? C*β/(a-γ)*((1-exp(-γ*t))/γ - (1-exp(-a*t))/a) : C*β/a*((1-exp(-γ*τ))/γ + exp(-a*t)*(1-exp(-(a-γ)τ))/(a-γ))
# where a = β + γ

C, γ, β, τ = [2., 0.1, 0.5, 15.]
a = β + γ 
x_A(t) = C/a*(1-exp(-a*t))
x_I(t)= 0<=t<=τ ? C*β/(a-γ)*((1-exp(-γ*t))/γ - (1-exp(-a*t))/a) : C*β/a*((1-exp(-γ*τ))/γ + exp(-a*t)*(1-exp((a-γ)τ))/(a-γ))


# Create a 
using DelaySSAToolkit
rate1 = [C,γ,β,γ,γ]
reactant_stoch = [[],[1=>1],[1=>1],[2=>1],[3=>1]]
net_stoch = [[1=>1],[1=>-1],[1=>-1,2=>1,3=>1],[2=>-1],[3=>-1]]
mass_jump = MassActionJump(rate1, reactant_stoch, net_stoch; scale_rates =false)
jumpsets = JumpSet((),(),nothing,[mass_jump])
# DelaySSA.var_to_jumps_map(2,mass_jump)

delay_trigger_affect! = function (de_chan, rng)
   append!(de_chan[1], τ)
   append!(de_chan[2], τ)
end
delay_trigger = Dict(3=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>-1],2=>[3=>-1]) # 1 代表 delay channel idx , value 代表 对 第 2 个 species -1 

delay_affect1! = function (de_chan, rng)
    i = rand(rng, 1:length(de_chan[1]))
    deleteat!(de_chan[1],i)
end
delay_affect2! = function (de_chan, rng)
    i = rand(rng, 1:length(de_chan[2]))
    deleteat!(de_chan[2],i)
end
delay_interrupt = Dict(4=>delay_affect1!,5=>delay_affect2!) 
delaysets = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)

#  vars = [1,2]
#  var_to_jumps = DelaySSA.var_to_jumps_map(2,mass_jump)
#  dep_rxs = reduce(vcat,[var_to_jumps[vars[i]] for i in eachindex(vars)])


u0 = [0,0,0]
tf = 30.
saveat = .1
# de_chan0 = Dict(1 =>[]) # No X_I for the initial delay channel
de_chan0 = [[],[]]
p = 0.
tspan = (0.,tf)
# aggregatoralgo = DelayMNRM()
aggregatoralgo = DelayRejection()
# aggregatoralgo = DelayDirect()
save_positions = (false,false)
dprob = DiscreteProblem(u0, tspan, p)
# jprob = JumpProblem( )
djprob = DelayJumpProblem(dprob, aggregatoralgo, jumpsets,delaysets,de_chan0,save_positions = (false,false))

#Fine tunning
# aggregator = djprob.jump_callback.discrete_callbacks[1].initialize

# integrator = DiffEqBase.__init(djprob,SSAStepper())
# integrator.cb.affect!

sol =@time solve(djprob, SSAStepper(),seed=10, saveat =.1, save_delay_channel = false)


using Plots, DiffEqBase
plot(sol, label = ["X_A" "X_I1" "X_I2"])

Sample_size = Int(1e4)
ens_prob = EnsembleProblem(djprob)
ens =@time solve(ens_prob,SSAStepper(),EnsembleThreads(),trajectories = Sample_size, saveat = .1, save_delay_channel =false)

using StatsBase
mean_A(t) = mean([ens[s](t)[1] for s in 1:Sample_size])
mean_I1(t) = mean([ens[s](t)[2] for s in 1:Sample_size])
mean_I2(t) = mean([ens[s](t)[3] for s in 1:Sample_size])
using Plots
timestamps = 0:0.1:30

plot(timestamps,x_A.(timestamps),linewidth=3)
plot!(timestamps,mean_A.(timestamps),linewidth=3,line=:dash)

plot(timestamps,x_I.(timestamps),linewidth=3)
plot!(timestamps,mean_I1.(timestamps),linewidth=3,line=:dash)
plot!(timestamps,mean_I2.(timestamps),linewidth=3,line=:dash, legend =:topleft)

