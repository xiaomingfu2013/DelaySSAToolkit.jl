# A telegraph model of stochastic delay 

## Model
According to [1], a telegraph gene expression model with stochastic delay assumes that the gene can switch between active $G$ and inactive $G^*$ states, transcribes nascent mRNA (denoted by  $N$) while in the active state which subsequently is removed a time $\tau$ later. Here the delay $\tau$ can be a random variable. The reaction scheme is given by
```math
G^*\xrightarrow{k_{\text{on}}} G,\\
G\xrightarrow{k_{\text{off}}}G^*,\\
G\xrightarrow{\rho}G+N.
```
Note that $G\xrightarrow{\rho}G+N$ will trigger $N\Rightarrow \emptyset$ after a delay time $\tau$. We set $k_{\text{on}}=0.0282$, $k_{\text{off}}=0.609$ and $\rho=2.11$.

### Markovian part
We first define the model using Catalyst (see [this example](tutorials.md) for more details about the construction of a reaction network).
```julia
using DiffEqJump, Catalyst, DelaySSAToolkit
using Random, Distributions
rn = @reaction_network begin
    kon, Goff --> Gon
    koff, Gon --> Goff
    ρ, Gon --> Gon + N
end kon koff ρ
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)
```
Then we set the Initial value and define a `DiscreteProblem`.
```julia
u0 = [1,0,0] # Gon, Goff, N
tf = 2000.
tspan = (0,tf)
dprob = DiscreteProblem(u0, tspan)
```
### Non-Markovian part
Unlike other examples, the elongation time $\tau$ is a random variable sampled from two different LogNormal distributions. We assume $\tau\sim \text{LogNormal}(0,2)+120$ and $\tau\sim \text{LogNormal}(1,\sqrt{2})+120$. For instance, we take $\tau\sim \text{LogNormal}(1,\sqrt{2})+120$ and define 
```julia
delay_trigger_affect! = function (integrator, rng)
    τ=rand(LogNormal(1,sqrt(2)))+120
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(3=>delay_trigger_affect!)
delay_complete = Dict(1=>[3=>-1]) 
delay_interrupt = Dict() 
delayjumpset = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)
```
To see how to define the `DelayJumpSet`, we refers to [this example](tutorials.md).
Thus, we can define the problem
```julia
de_chan0 = [[]]
djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(false,false))
```
## Visualisation
We simulate $10^5$ trajectories and calculate the probability distribution.
```julia
ensprob = EnsembleProblem(djprob)
@time ens = solve(ensprob, SSAStepper(), EnsembleThreads(), trajectories=10^5)
```
![stochastic_delay1](../assets/stochastic_delay1.svg)

If we change the stochastic delay to $\tau\sim \text{LogNormal}(0,2)+120$, we can see the zero-inflated mode disappeared.
![stochastic_delay2](../assets/stochastic_delay2.svg)
## References

[1] Qingchao Jiang, Xiaoming Fu, Shifu Yan, Runlai Li, Wenli Du, Zhixing Cao, Feng Qian, Ramon Grima, "Neural network aided approximation and parameter inference of non-Markovian models of gene expression". Nature communications, (2021) 12(1), 1-12. [https://doi.org/10.1038/s41467-021-22919-1](https://doi.org/10.1038/s41467-021-22919-1)