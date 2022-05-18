using DelaySSAToolkit
using Catalyst
using Statistics
using Test


reltol = 1e-1
@parameters t
@variables X(t)
burst_sup = 30
a, b = [0.0282, 3.46]
rxs = [Reaction(a*b^i/(1+b)^(i+1),nothing,[X],nothing,[i]) for i in 1:burst_sup]
rxs = vcat(rxs)
@named rs = ReactionSystem(rxs,t,[X],[])

# convert the ReactionSystem to a JumpSystem
jumpsys = convert(JumpSystem, rs, combinatoric_ratelaws=false)

u0 = [0]
de_chan0 = [[]]
tf = 200.
tspan = (0,tf)
dprob = DiscreteProblem(jumpsys,u0,tspan)
τ = 130.

delay_trigger = Dict([Pair(i, [1=>fill(τ, i)]) for i in 1:burst_sup])
delay_complete = Dict(1=>[1=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)



algos = [DelayRejection(), DelayDirect(),DelayMNRM(),  DelayDirectCR()]
timestamps = [10, 20, 50, 200]

bursty_mean(t) = a*b*min(t,τ)
bursty_var(t) = 2*a*b^2*min(t,τ) + a*b*min(t,τ)

# Check mean and variance
samplesize = Int64(5e4)
@testset for alg in algos
    jprob = DelayJumpProblem(jumpsys, dprob, alg, delayjumpset, de_chan0, save_positions=(false, false))
    ensprob = EnsembleProblem(jprob)
    @info "Testing method $(alg)"
    @time ens = solve(ensprob, SSAStepper(), EnsembleSerial(), trajectories=samplesize, saveat = timestamps)
    for idx in eachindex(timestamps)
    sol_end = reduce(vcat, [ens[i].u[idx+1] for i in 1:samplesize])
    t_ = timestamps[idx]
    @test mean(sol_end) ≈ bursty_mean(t_) atol = reltol*bursty_mean(t_) 
    @test var(sol_end) ≈ bursty_var(t_) atol = reltol*bursty_var(t_)
    end
end

