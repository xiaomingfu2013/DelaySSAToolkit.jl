# The model is 
# C: 0 -> X_A, 
# γ: X_A -> 0 
# β: X_A ->X_I : trigger X_I  -> 0 after τ time
# γ: X_I ->0 : delay_interrupt

# According to Lafuerza, Toral 2011 the exact solution is obtainable
# if we denote x_A(t) to be the mean value of X_A at time t, and x_I(t) the mean value of X_I at time t, then 
# x_A(t)= C/a (1-e^(-at))  
# x_I(t)= 0<=t<=τ ? C*β/(a-γ)*((1-exp(-γ*t))/γ - (1-exp(-a*t))/a) : C*β/a*((1-exp(-γ*τ))/γ + exp(-a*t)*(1-exp(-(a-γ)τ))/(a-γ))
# where a = β + γ


using DiffEqJump, Catalyst
using DelaySSAToolkit 

rn = @reaction_network begin
   C, 0 --> Xₐ
   γ, Xₐ --> 0
   β, Xₐ --> Xᵢ 
   γ, Xᵢ --> 0
end C γ β

jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)

u0 = [0, 0]
tf = 30.
saveat = .1
de_chan0 = [[]]
C, γ, β = [2., 0.1, 0.5]
p = [C, γ, β]
tspan = (0.,tf)
aggregatoralgo = DelayRejection()
# aggregatoralgo = DelayMNRM()
# aggregatoralgo = DelayDirect()
# aggregatoralgo = DelayDirectCR()
dprob = DiscreteProblem(u0, tspan, p)
djprob = DelayJumpProblem(jumpsys, dprob, aggregatoralgo,  delaysets, de_chan0, save_positions = (false, false))


τ = 15.
delay_trigger_affect! = function (integrator, rng)
   append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(3=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>-1]) 
delay_affect! = function (integrator, rng)
    i = rand(rng, 1:length(integrator.de_chan[1]))
    deleteat!(integrator.de_chan[1],i)
end
delay_interrupt = Dict(4=>delay_affect!) 
delaysets = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)

sol =@time solve(djprob, SSAStepper(), seed = 2, saveat =.1, save_delay_channel = true)
sol =@time solve(djprob, SSAStepper(), seed = 2, save_delay_channel = true)

using Plots, DiffEqBase; theme(:vibrant)
ens_prob = EnsembleProblem(djprob)
Sample_size = Int(10^4)
@time ens = solve(ens_prob, SSAStepper(),EnsembleThreads(),trajectories = Sample_size, saveat = .1)
plot(ens[1], label = ["X_A" "X_I"], fmt =:svg)
# savefig("docs/src/assets/delay_degradation1.svg")


a = β + γ 
x_A(t) = C/a*(1-exp(-a*t))
x_I(t)= 0<=t<=τ ? C*β/(a-γ)*((1-exp(-γ*t))/γ - (1-exp(-a*t))/a) : C*β/a*((1-exp(-γ*τ))/γ + exp(-a*t)*(1-exp((a-γ)τ))/(a-γ))


using StatsBase
mean_A(t) = mean([ens[s](t)[1] for s in 1:Sample_size])
mean_I(t) = mean([ens[s](t)[2] for s in 1:Sample_size])

timestamps = 0:0.1:tf
plot(timestamps,x_A.(timestamps),linewidth=3, label = "X_A Exact", xlabel = "Time", ylabel = "# of X_A")
plot!(timestamps,mean_A.(timestamps),linewidth=3,line=:dash, label =  "X_A SSA")
plot!(timestamps,x_I.(timestamps),linewidth=3, label = "X_I Exact")
plot!(timestamps,mean_I.(timestamps),linewidth=3,line=:dash, legend = :topleft, label = "X_I SSA")
# savefig("docs/src/assets/delay_degradation2.svg")


x_A.(timestamps)
mean_A.(timestamps)
x_I.(timestamps)
mean_I.(timestamps)
