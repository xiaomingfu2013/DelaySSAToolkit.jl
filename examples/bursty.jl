## This example shows an application of `DelaySSA.jl` for a bursty model with delay

## A Bursty model with delay is described as 
# ab^n/(1+b)^(n+1): 0 -> n P, which triggers n P to degrade after delay time τ

using DelaySSAToolkit
using Catalyst, DiffEqJump


begin # construct reaction network
    @parameters a b t
    @variables X(t)
    burst_sup = 30
    rxs = [Reaction(a*b^i/(1+b)^(i+1),nothing,[X],nothing,[i]) for i in 1:burst_sup]
    rxs = vcat(rxs)
    @named rs = ReactionSystem(rxs,t,[X],[a,b])
end

# convert the ReactionSystem to a JumpSystem
jumpsys = convert(JumpSystem, rs, combinatoric_ratelaws=false)

u0 = [0]
de_chan0 = [[]]
tf = 200.
tspan = (0,tf)
timestamp = 0:1:tf
ps = [0.0282, 3.46]
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
τ = 130.

delay_trigger_affect! = []
for i in 1:burst_sup
    push!(delay_trigger_affect!, function (de_chan, rng)
    append!(de_chan[1], fill(τ, i))
    end)
end
delay_trigger_affect!
delay_trigger = Dict([Pair(i, delay_trigger_affect![i]) for i in 1:burst_sup])
delay_complete = Dict(1=>[1=>-1])
delay_interrupt = Dict()

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
jprob = DelayJumpProblem(jumpsys, dprob, DelayDirect(), delayjumpset, de_chan0, save_positions=(true,true))
sol = solve(jprob, SSAStepper(), seed= 1)

using Plots, DiffEqBase
ensprob = EnsembleProblem(jprob)
@time ens = solve(ensprob, SSAStepper(), EnsembleThreads(),saveat=timestamp, trajectories=10^4)

# Check with the exact probability distribution
using TaylorSeries
function taylor_coefficients(NT::Int,at_x,gen::Function)
    Q = zeros(NT)
    taylor_gen = taylor_expand(u->gen(u),at_x,order=NT)
    for j in 1:NT
        Q[j] = taylor_gen[j-1]
    end
    return Q
end
function delay_bursty(params,NT::Int)
    a, b, τ, t = params
    gen1(u) = exp(a*b*min(τ,t)*u/(1-b*u))
    taylor_coefficients(NT,-1,gen1)
end


using DifferentialEquations.EnsembleAnalysis
theme(:vibrant)
sol_end = componentwise_vectors_timepoint(ens, tf)
histogram(sol_end,bins=0:1:60,normalize=:pdf, label = "SSA",fillalpha = 0.6, linecolor = :orange)
sol_exact = delay_bursty([ps;130;200], 60)
plot!(0:59,sol_exact, linewidth = 3, label = "Exact solution")