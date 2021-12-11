using Catalyst
using Distributions, Random
using DiffEqJump
using DelaySSAToolkit

rn = @reaction_network begin
    A, 0 --> X
    B, Y --> 0
end A B

u0 = [0,0]
de_chan0 = [[]]
tf = 100.
tspan = (0,tf)
ps = [10., 1.] # dummy variables, later on will be drawn from a Gamma distribution 
τ = 1 # dummy variables, later on will be drawn from a Gamma distribution 
delay_trigger = Dict(1=>[1=>τ])
delay_complete = Dict(1=>[1=>-1, 2=>1])
delay_interrupt = Dict()


jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
# use ensemble problem 
algo = DelayRejection()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
djprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0, save_positions=(false,false))

function prob_func(prob, i ,repeat)
    rng = Random.seed!(i)
    α = rand(rng, Gamma(63,1/9)) # the Gamma distribution uses shape α, scale θ, Gamma(α,θ). In the paper, Gamma distribution uses shape α and rate  β. One needs to set the inverse.
    β = rand(rng, Gamma(10.,1/10.))
    A = rand(rng, Gamma(8,1/0.23))
    B = rand(rng, Gamma(9,1/625))
    τ  = rand(rng, Gamma(α,1/β))
    delay_trigger_i = Dict(1=>[1=>τ])
    prob.delayjumpsets.delay_trigger = delay_trigger_i
    @. prob.prob.p = [A, B]
    DiffEqJump.update_parameters!(prob.massaction_jump, [A, B])
    prob
end

function prob_func_simple(prob, i, repeat)
    rng = Random.seed!(i)
    α = rand(rng, Gamma(63,1/9)) # the Gamma distribution uses shape α, scale θ, Gamma(α,θ). In the paper, Gamma distribution uses shape α and rate  β. One needs to set the inverse.
    β = rand(rng, Gamma(10.,1/10.))
    
    A = rand(rng, Gamma(8,1/0.23))
    B = rand(rng, Gamma(9,1/625))
    τ  = rand(rng, Gamma(α,1/β))
    remake(prob, p = [A,B], delay_trigger=Dict(1=>[1=>τ]))
end
prob1 = EnsembleProblem(djprob, prob_func = prob_func)
prob2 = EnsembleProblem(djprob, prob_func = prob_func_simple)

@time ens = solve(prob1, SSAStepper(), EnsembleThreads(), trajectories = 10^3, saveat= 1.)
@time ens2 = solve(prob2, SSAStepper(), EnsembleThreads(), trajectories = 10^3, saveat=1.)

using Plots
plot(ens[1])
plot!(ens2[1])

using DifferentialEquations.EnsembleAnalysis
last_slice = componentwise_vectors_timepoint(ens, tf)
histogram(last_slice[2], bins=0:1:30, normalize=:pdf)

delayjumpset


new_djprob = remake(djprob, delay_trigger = Dict(2=>[1=>2]), p = [1,2.], de_chan0 = [[1.]], u0 = [1,0], delay_complete = Dict(2=>[2=>-1]), delay_interrupt_set = [2])

new_djprob.prob.p
new_djprob.massaction_jump
new_djprob.de_chan0
new_djprob.prob.u0
new_djprob.delayjumpsets

djprob.delayjumpsets
