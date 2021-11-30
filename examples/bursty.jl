## This example shows an application of `DelaySSA.jl` for a bursty model with delay


## A Bursty model with delay is described as 
# ab^n/(1+b)^(n+1): 0 -> n P, which triggers n P to degrade after delay time τ

using DelaySSA
using Catalyst
using Plots, DiffEqBase

begin # construct reaction network
    @parameters a b t
    @variables X(t)
    burst_sup = 30
    rxs = [Reaction(a*b^i/(1+b)^(i+1),nothing,[X],nothing,[i]) for i in 1:burst_sup]
    rxs = vcat(rxs)
    @named rs_new = ReactionSystem(rxs,t,[X],[a,b])
end;
rs_new.eqs
jumpsys = convert(JumpSystem, rs_new, combinatoric_ratelaws=false)

u0 = [0]
de_chan0 = [[]]
tf = 20.
tspan = (0,tf)
timestamp = 0:1:tf
ps = [0.0282, 3.46]
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
# append!([],fill(1,10))
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

# jprob = DelayJumpProblem(jumpsys, dprob, DelayDirect(), delayjumpset, de_chan0, save_positions=(false,false))
jprob.massaction_jump

equations(jumpsys)