using JumpProcesses, DelaySSAToolkit
using Random, Distributions

rn = @reaction_network begin
    kon, Goff --> Gon
    koff, Gon --> Goff
    ρ, Gon --> Gon + N
end
jumpsys = convert(JumpSystem, rn; combinatoric_ratelaws=false)

u0 = [1, 0, 0]
de_chan0 = [[]]
tf = 2000.0
tspan = (0, tf)
p = [0.0282, 0.609, 2.11]
dprob = DiscreteProblem(u0, tspan, p)

delay_trigger_affect! = function (integrator, rng)
    τ = rand(LogNormal(1, sqrt(2))) + 120
    return append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(3 => delay_trigger_affect!)
delay_complete = Dict(1 => [3 => -1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

djprob = DelayJumpProblem(
    jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0; save_positions=(false, false)
)

ensprob = EnsembleProblem(djprob)
@time ens = solve(ensprob, SSAStepper(), EnsembleThreads(); trajectories=10^5)

using DifferentialEquations.EnsembleAnalysis, Plots;
theme(:vibrant);
last_slice = componentwise_vectors_timepoint(ens, tf)
histogram(
    last_slice[3];
    bins=0:1:60,
    normalize=:pdf,
    label="LogNormal(1,sqrt(2))",
    fillalpha=0.6,
    linecolor=:orange,
    fmt=:svg,
    xlabel="# of products",
    ylabel="Probability",
)
# savefig("docs/src/assets/stochastic_delay1.svg")

delay_trigger_affect2! = function (integrator, rng)
    τ = rand(LogNormal(0, 2)) + 120
    return append!(integrator.de_chan[1], τ)
end
delay_trigger2 = Dict(3 => delay_trigger_affect2!)
delayjumpset2 = DelayJumpSet(delay_trigger2, delay_complete, delay_interrupt)

djprob2 = DelayJumpProblem(
    jumpsys, dprob, DelayRejection(), delayjumpset2, de_chan0; save_positions=(false, false)
)

@time ens2 = solve(
    EnsembleProblem(djprob2), SSAStepper(), EnsembleThreads(); trajectories=10^5
)

last_slice2 = componentwise_vectors_timepoint(ens2, tf)
histogram(
    last_slice2[3];
    bins=0:1:60,
    normalize=:pdf,
    label="LogNormal(0,2)",
    fillalpha=0.6,
    linecolor=:orange,
    fmt=:svg,
    xlabel="# of products",
    ylabel="Probability",
)
# savefig("docs/src/assets/stochastic_delay2.svg")
