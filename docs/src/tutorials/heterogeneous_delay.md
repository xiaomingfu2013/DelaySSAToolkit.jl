# A stochastic birth-death model with heterogeneous fixed delays

## Model definition
We study the following model
```math
A_n: \emptyset \rightarrow X \text{ triggers } X \Rightarrow Y \text{ after $\tau$ time}\\
B_n: Y \rightarrow \emptyset
```
This model ....

### Markovian part
```julia
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
ps = [10., 1.] # A, B. The initial conditions dummy variables, later on will be drawn from a Gamma distribution 
```

### Non-Markovian part

```julia
τ = 1 # dummy variable, later on will be drawn from a Gamma distribution 
delay_trigger = Dict(1=>[1=>τ])
delay_complete = Dict(1=>[1=>-1, 2=>1])
delay_interrupt = Dict()
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
jprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(false,false))
```

## Define `EnsembleProblem` 
```julia
function prob_func(prob, i ,repeat)
    rng = Random.seed!(i)
    α = rand(rng, Gamma(63,1/9)) # the Gamma distribution uses shape α, scale θ, Gamma(α,θ). In the paper, Gamma distribution uses shape α and rate  β. One needs to set the inverse.
    β = rand(rng, Gamma(10.,1/10.))
    A = rand(rng, Gamma(8,1/0.23))
    B = rand(rng, Gamma(9,1/625))
    τ = rand(rng, Gamma(α,1/β))
    delay_trigger_i = Dict(1=>[1=>τ])
    prob.delayjumpsets.delay_trigger = delay_trigger_i
    @. prob.prob.p = [A, B]
    prob
end
ensprob = EnsembleProblem(jprob, prob_func = prob_func)
@time ens = solve(ensprob, SSAStepper(), EnsembleThreads(), trajectories = 10^4, saveat = 1.)
```

For distributed delay, one can easily use the following setting.
## Define `EnsembleProblem` 
```julia
function prob_func(prob, i ,repeat)
    rng = Random.seed!(i)
    α = rand(rng, Gamma(63,1/9)) # the Gamma distribution uses shape α, scale θ, Gamma(α,θ). In the paper, Gamma distribution uses shape α and rate  β. One needs to set the inverse.
    β = rand(rng, Gamma(10.,1/10.))
    A = rand(rng, Gamma(8,1/0.23))
    B = rand(rng, Gamma(9,1/625))
    τ = rand(rng, Gamma(α,1/β))
    delay_trigger_i = Dict(1=>[1=>τ])
    prob.delayjumpsets.delay_trigger = delay_trigger_i
    @. prob.prob.p = [A, B]
    prob
end
ensprob = EnsembleProblem(jprob, prob_func = prob_func)
```

## Visualisation
```julia
@time ens = solve(ensprob, SSAStepper(), EnsembleThreads(), trajectories = 10^4, saveat = 1.)
```
```julia-repl
30.610441 seconds (196.13 M allocations: 39.616 GiB, 22.67% gc time)
EnsembleSolution Solution of length 10000 with uType
```
