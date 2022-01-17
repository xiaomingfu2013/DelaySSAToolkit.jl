using DelaySSAToolkit, Catalyst
using LaTeXStrings
using Plots
theme(:vibrant)

rn = @reaction_network begin
    C, 0 --> Xₐ
    γ, Xₐ --> 0
    β, Xₐ --> Xᵢ₁ + Xᵢ₂
    γ, Xᵢ₁ --> 0
    γ, Xᵢ₂ --> 0
 end C γ β
 jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)
 
C, γ, β, τ = [2., 0.1, 0.5, 15.]
p=[C, γ, β]
u0 = [0, 0, 0]
tf = 30.
saveat = .1
de_chan0 = [[],[]]
tspan = (0.,tf)
dprob = DiscreteProblem(u0, tspan,p)

# delay_trigger_affect! = function (integrator, rng)
#    append!(integrator.de_chan[1], τ)
#    append!(integrator.de_chan[2], τ)
# end
# delay_trigger = Dict(3=>delay_trigger_affect!)
delay_trigger = Dict(3=>[1=>τ,2=>τ])
delay_complete = Dict(1=>[2=>-1],2=>[3=>-1]) 

delay_affect1! = function (integrator, rng)
    i = rand(rng, 1:length(integrator.de_chan[1]))
    deleteat!(integrator.de_chan[1],i)
end
delay_affect2! = function (integrator, rng)
    i = rand(rng, 1:length(integrator.de_chan[2]))
    deleteat!(integrator.de_chan[2],i)
end
delay_interrupt = Dict(4=>delay_affect1!,5=>delay_affect2!) 
delayjumpset = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)

djprob = DelayJumpProblem(jumpsys,dprob, DelayRejection(),  delayjumpset, de_chan0, save_positions=(false,false))
sol =@time solve(djprob, SSAStepper(),seed=10, saveat =.1, save_delay_channel = false)
plot(sol, label = [L"X_A" L"X_{I1}" L"X_{I2}"],legend =:topleft,xlabel="time",ylabel="Molecular Number")

timestamps = 0:0.1:tf
a = β + γ 
x_A(t) = C/a*(1-exp(-a*t))
x_I(t)= 0<=t<=τ ? C*β/(a-γ)*((1-exp(-γ*t))/γ - (1-exp(-a*t))/a) : C*β/a*((1-exp(-γ*τ))/γ + exp(-a*t)*(1-exp((a-γ)τ))/(a-γ))
using Plots
plot(timestamps,x_A.(timestamps),linewidth=3,label=L"X_A",xlabel="time",ylabel="Mean Value")
plot!(timestamps,x_I.(timestamps),linewidth=3,label=L"X_I")

using StatsBase, DifferentialEquations.EnsembleAnalysis
Sample_size = Int(1e4)
ens_prob = EnsembleProblem(djprob)
ens =@time solve(ens_prob, SSAStepper(),EnsembleThreads(),trajectories = Sample_size, saveat = .1)

tsm(t) = timepoint_mean(ens,t)

mean_A(t) = tsm(t)[1]
mean_I1(t) = tsm(t)[2]
mean_I2(t) = tsm(t)[3]


plot(timestamps,x_A.(timestamps),linewidth=3,label=L"$X_A$ Exact",xlabel="time",ylabel="Mean Value")
plot!(timestamps,mean_A.(timestamps),linewidth=3,line=:dash,label=L"$X_A$ Delay SSA")
plot!(timestamps,x_I.(timestamps),linewidth=3,label=L"$X_I$ Exact")
plot!(timestamps,mean_I1.(timestamps),linewidth=3,line=:dash,label=L"$X_{I1}$ Delay SSA")
plot!(timestamps,mean_I2.(timestamps),linewidth=3,line=:dash, legend =:topleft,label=L"$X_{I2}$ Delay SSA")


