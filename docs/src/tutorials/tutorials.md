# [Tutorials](@id seir_model)

This tutorial is aims to explain how to use DelaySSAToolkit to define chemical reaction models, solve the problem and visualize the results. To demonstrate the functionalities, we will consider a specific case from epideimc modelling as follows
```math
S+I\xrightarrow{\rho}E+I\\
I\stackrel{r}{\rightarrow}R
```
and $S+I\xrightarrow{\rho} E+I$ will trigger $E\Rightarrow I$ after $\tau$ time, where $S$, $I$ and $R$ are the susceptible, infected and removed populations. $E$ represents the exposed population. It means, wtih rate $\rho$, a susceptible contacted by an infected will become an individual that is exposed to the disease and then it takes certain amount of time delay $\tau$ to become an infected inidividual.

# Model

What differs from the Markov process that can be modelled via SSA is the introduction of **delay reactions**. To show how we incorporate the delay reactions into the Markovian system, we first need to define the Markovian part and then its non-Markovian part. These two parts mainly form a `DelayJumpProblem`. Here we show two routes to define our delay system, one way is based on `JumpSystem`, `DiscreteProblem` and `DelayJumpSet`, the other is based on `JumpSet`, `DiscreteProblem` and `DelayJumpSet`.

## First route: `JumpSystem + DiscreteProblem + DelayJumpSet`
### [Markovian part](@id Markovian_part)
[Catalyst.jl](https://github.com/SciML/Catalyst.jl) provides a comprehensive interface to modelling reaction networks in Julia and can be used to construct models fully-compatible with DelaySSAToolkit. For more details on how to construct a reaction network, we recommend reading [Catalyst's tutorial](https://catalyst.sciml.ai/stable/tutorials/using_catalyst/). In our example, the model can be defined as:
```julia
rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r
```
We can easily obtain `Jumpsystem` from the reaction network `rn` that is previously defined using Catalys interface.

```julia
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
```
where `combinatoric_ratelaws` is an optional parameter that specifies whether the rate constants correspond to stochastic rate constants in the sense used by Gillespie, and hence need to be rescaled. The default, `combinatoric_ratelaws=true`, corresponds to rescaling the passed in rate constants. The default behavior is to assume rate constants correspond to stochastic rate constants in the sense used by Gillespie [1]. This means that for a reaction such as $2A \overset{k}{\rightarrow} B$, the jump rate function constructed by `MassActionJump` would be `k*A*(A-1)/2!`. For a trimolecular reaction like $3A \overset{k}{\rightarrow} B$ the rate function would be `k*A*(A-1)*(A-2)/3!`. To *avoid* having the reaction rates rescaled (by `1/2` and `1/6` for these two examples), one can pass the `JumpSystem` constructor the optional named parameter `combinatoric_ratelaws=false` see [Reaction rate laws used in simulations](https://catalyst.sciml.ai/stable/tutorials/using_catalyst/#Reaction-rate-laws-used-in-simulations) for details.

With the initial conditions, we can then define `DiscreteProblem`
```julia
u0 = [999,1,0,0] # S, I, E, R
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2] # parameters for ρ, r
τ = 20.
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
```
where `DiscreteProblem` inputs `jumpsys`, and the initial condition of reactants `u0` , the simulation timespan `tspan` and the reaction rates `ps`.

### Non-Markovian part
The non-Markovian part consists of three elements:
- delay trigger reactions: those reactions in the [Markovian part](@ref Markovian_part) that trigger the change of the state of the delay channels or/and the state of the reactants upon initiation.
- delay interrupt reactions: those reactions in the [Markovian part](@ref Markovian_part) that change the state of the delay channels or/and the state of the reactants in the middle of on-going delay reactions.
- delay complete reactions: those reactions that are initiated by delay trigger reactions and change the state of the delay channels or/and the state of the reactants upon completion.
  
With these three definitions in mind and based on this particular example, we define the `DelayJumpSet` by
```julia
delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_interrupt = Dict()
delay_complete = Dict(1=>[2=>1, 3=>-1])
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
```
- `delay_trigger::Dict`  A dictionary that contains
  - Keys: Indices of reactions defined in `JumpSystem` that can trigger the delay reaction see [indices in the `JumpSystem`](@ref indice_notice). Here we have the first reaction $S+I\Rightarrow E+ I$ that will trigger the transfer from $E$ to $I$ after time $\tau$, hence the key here is `1`.
  
  - Values: A update function that determines how to update the delay channel. In this example, once the delay reaction is trigged, the first delay channel will be added a delay time $\tau$.
  
- `delay_interrupt::Dict`
  - There are no delay interrupt reactions in this example so we set `delay_interrupt = Dict()`.
- ```delay_complete::Dict``` 
  - Keys: Indices of delay channels.
  - Values: A vector of `Pair`s, mapping species index to net change of stoichiometric coefficient.

We refer to [Defining a `DelayJumpSet`](bursty.md/#Defining-a-DelayJumpSet) and [Defining a `DelayJumpSet`(birth-death example)](delay_degradation.md/#Defining-a-DelayJumpSet) for more details.

### [Remark on Reaction Indices](@id indice_notice) 
!!! warning
    `JumpSystem` might change the order of the reactions that is arranged in your reaction network. Internally, all MassActionJumps are ordered before ConstantRateJumps (with the latter internally ordered in the same order they were passed in). The same principle applies for the construction of `JumpSet`.

At last, we can define the `DelayJumpProblem` by 
```julia
jprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(true,true))
```
where `DelayJumpProblem` inputs `jumpsys`,`DiscreteProblem`, `DelayJumpSet`, the algorithm we choose and the initial condition of the delay channel `de_chan0`.

Let's put aside our `DelayJumpProblem` for the moment and see how to define the delay system in another way. 

## Second route: `JumpSet + DiscreteProblem + DelayJumpSet`
Now we explain how to define the `DelayJumpProblem` in another way. To that aim, we should first define the parameters and the mass-action jump (see [Defining a Mass Action Jump](https://diffeq.sciml.ai/stable/types/jump_types/#Defining-a-Mass-Action-Jump) for details) and construct `Jumpset`.
### Markovian part
```julia 
ρ, r = [1e-4, 1e-2]
rate1 = [ρ, r]
reactant_stoich = [[1=>1,2=>1],[2=>1]]
net_stoich = [[1=>-1,3=>1],[2=>-1,4=>1]]
mass_jump = MassActionJump(rate1, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,mass_jump)
```
The `JumpSet` consists of four inputs, namely variable jumps, constant rate jumps, regular jumps and mass-action jumps. As far as discrete stochastic simulation is concerned, we only focus on constant rate jumps and mass-action jumps which is the second and fourth inputs of `JumpSet` (see [different jump types](https://diffeq.sciml.ai/stable/types/jump_types/) for more details). Here we only have two mass-action jumps that are wrapped in `mass_jump`.

Then we initialise the problem by setting
```julia
u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
τ = 20.
```
As before, we can define the `DiscreteProblem`
```Julia
dprob = DiscreteProblem(u0, tspan)
```
### Non-Markovian part
In the same way, we can define the  `DelayJumpSet` by
```julia
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
```
Now we can define the problem
```julia 
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions=(true,true)).
```
At last, we can solve the problem and visualize it
```julia
sol = solve(djprob, SSAStepper())
```
![seir](../assets/seir.svg)
