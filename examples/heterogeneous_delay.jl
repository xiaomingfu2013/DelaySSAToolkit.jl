using Catalyst
using Distributions, Random
using DelaySSAToolkit

rn = @reaction_network begin
    A, 0 --> X
    B, Y --> 0
end A B

u0 = [0, 0]
de_chan0 = [[]]
tf = 100.0
tspan = (0, tf)
ps = [10.0, 1.0] # dummy variables, later on will be drawn from a Gamma distribution 
τ = 1 # dummy variables, later on will be drawn from a Gamma distribution 
delay_trigger = Dict(1 => [1 => τ])
delay_complete = Dict(1 => [1 => -1, 2 => 1])
delay_interrupt = Dict()

jumpsys = convert(JumpSystem, rn; combinatoric_ratelaws=false)
dprob = DiscreteProblem(jumpsys, u0, tspan, ps)
# use ensemble problem 
algo = DelayRejection()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
djprob = DelayJumpProblem(
    jumpsys, dprob, algo, delayjumpset, de_chan0; save_positions=(false, false)
)

function prob_func1(prob, i, repeat)
    rng = Random.seed!(i)
    A = rand(rng, Gamma(8, 1 / 0.23))
    B = rand(rng, Gamma(9, 1 / 625))
    τ = rand(rng, Gamma(7, 1))
    return remake(prob; p=[A, B], delay_trigger=Dict(1 => [1 => τ]))
end
function prob_func2(prob, i, repeat)
    rng = Random.seed!(i)
    α = rand(rng, Gamma(63, 1 / 9)) # the Gamma distribution uses shape α, scale θ, Gamma(α,θ). In the paper, Gamma distribution uses shape α and rate  β. One needs to set the inverse.
    β = rand(rng, Gamma(10.0, 1 / 10.0))
    A = rand(rng, Gamma(8, 1 / 0.23))
    B = rand(rng, Gamma(9, 1 / 625))
    τ = rand(rng, Gamma(α, 1 / β))
    return remake(prob; p=[A, B], delay_trigger=Dict(1 => [1 => τ]))
end

prob1 = EnsembleProblem(djprob; prob_func=prob_func1)
@time ens1 = solve(prob1, SSAStepper(), EnsembleThreads(); trajectories=40, saveat=1.0)
Y_evolution = [ens1[i][2, :] for i in 1:40]

using DifferentialEquations.EnsembleAnalysis
using Plots;
theme(:vibrant);
plot(
    Y_evolution;
    linealpha=0.4,
    legend=false,
    linewidth=2,
    fmt=:svg,
    xlabel="Time",
    ylabel="# of Y individuals",
)
savefig("docs/src/assets/heterogeneous_delay1.svg")

prob2 = EnsembleProblem(djprob; prob_func=prob_func2)
@time ens2 = solve(prob2, SSAStepper(), EnsembleThreads(); trajectories=10^4)

last_slice = componentwise_vectors_timepoint(ens2, tf)
histogram(
    last_slice[1];
    bins=0:20:1000,
    normalize=:pdf,
    xlabel="# of X individuals",
    ylabel="Probability",
    fillalpha=0.6,
    linecolor=:orange,
    legend=false,
)
# savefig("docs/src/assets/heterogeneous_delay2.svg")
