# A birth-death example with delay degradation

## Model definition

The model is defined as follows: 1. $C:\emptyset \rightarrow X_A$; 2. $\gamma : X_A \rightarrow \emptyset$; 3. $\beta : X_A \rightarrow  X_I$, which triggers $X_I\Rightarrow \emptyset$ after $\tau$ time; 4. $\gamma: X_I \rightarrow \emptyset$, which causes the delay channel to change its state during a schduled delay reaction.

This example is studied by Lafuerza and Toral in [1], where one can solve the solution analytically. If we denote $\langle X_A\rangle(t)$ to be the mean value of $X_A$ at time $t$, and $\langle X_I\rangle(t)$ the mean value of $X_I$ at time $t$, then
```math
\langle X_A\rangle(t)= \frac{C}{a}( 1-e^{-at} ),\quad \langle X_I\rangle(t) = \begin{cases}
\frac{C\beta}{a-γ}\big[\frac{1-e^{-γt}}{γ}-\frac{1-e^{-at}}{a}\big]，& t \in [0,\tau]\\
\frac{C\beta}{a}\Big[\frac{1-e^{-γτ}}{γ}+\frac{(1-e^{\tau(a-γ)})}{a-γ}e^{-at}\Big], & t \in (\tau,\infty)
\end{cases}
```
where $a = β + γ$.

We first define the parameters and the mass-action jump (see [Defining a Mass Action Jump](https://diffeq.sciml.ai/stable/types/jump_types/#Defining-a-Mass-Action-Jump) for details)

```julia
C, γ, β, τ = [2., 0.1, 0.5, 15.]
rate1 = [C,γ,β,γ]
reactant_stoich = [[],[1=>1],[1=>1],[2=>1]]
net_stoich = [[1=>1],[1=>-1],[1=>-1,2=>1],[2=>-1]]
mass_jump = MassActionJump(rate1, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,[mass_jump])
```

- `rates`  A vector of rates of reactions.
- `reactant_stoch` is a vector whose `k`th entry is the reactant stoichiometry of the `k`th reaction. The reactant stoichiometry for an individual reaction is assumed to be represented as a vector of `Pair`s, mapping species id to stoichiometric coefficient.
- `net_stoch`  is assumed to have the same type as `reactant_stoich`; a vector whose `k`th entry is the net stoichiometry of the `k`th reaction. The net stoichiometry for an individual reaction is again represented as a vector of `Pair`s, mapping species id to the net change in the species when the reaction occurs.
- `scale_rates` is an optional parameter that specifies whether the rate constants correspond to stochastic rate constants in the sense used by Gillespie, and hence need to be rescaled. *The default, `scale_rates=true`, corresponds to rescaling the passed in rate constants.* When using `MassActionJump` the default behavior is to assume rate constants correspond to stochastic rate constants in the sense used by Gillespie (J. Comp. Phys., 1976, 22 (4)). This means that for a reaction such as $2A \overset{k}{\rightarrow} B$, the jump rate function constructed by `MassActionJump` would be `k*A*(A-1)/2!`. For a trimolecular reaction like $3A \overset{k}{\rightarrow} B$ the rate function would be `k*A*(A-1)*(A-2)/3!`. To *avoid* having the reaction rates rescaled (by `1/2` and `1/6` for these two examples), one can pass the `MassActionJump` constructor the optional named parameter `scale_rates=false`
- `mass_jump`  Define mass-action jumps
- `jumpsets`  Wrap up the reactions into one jumpset.

### Defining a `DelayJumpSet`

Then we turn to the definition of delay reactions

```julia
delay_trigger_affect! = function (de_chan, rng)
   append!(de_chan[1], τ)
end
delay_trigger = Dict(3=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>-1]) 

delay_affect! = function (de_chan, rng)
    i = rand(rng, 1:length(de_chan[1]))
    deleteat!(de_chan[1],i)
end
delay_interrupt = Dict(4=>delay_affect!) 
delaysets = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)
```

- `delay_trigger`  
  - Keys: Indices of reactions defined in `jumpset` that can trigger the delay reaction. Here we have the 3rd reaction $\beta: X_A \rightarrow X_I$ that will trigger the $X_I$ to degrade after time $\tau$.
  - Values: A update function that determines how to update the delay channel. In this example, once the delay reaction is trigged, the delay channel 1 (which is the channel for $X_I$) will be added a delay time $\tau$.			
- `delay_interrupt`
  - Keys: Indices of reactions defined in `jumpset` that can cause the change in delay channel. In this example, the 4th reaction $\gamma : X_I \rightarrow \emptyset$ will change the schduled delay reaction to change its state immediately.
  - Values: A update function that determines how to update the delay channel. In this example, once a delay-interrupt reaction happens, any of the reactants $X_I$ that is supposed to leave the system after time $\tau$ can be degraded immediately.  
- `delay_complete` 
  - Keys: Indices of delay channel. Here the 1st delay channel corresponds to $X_I$.
  - Values: A vector of `Pair`s, mapping species id to net change of stoichiometric coefficient.

Now we can initialise the problem by setting 

```julia
u0 = [0, 0]
tf = 30.
saveat = .1
de_chan0 = [[]]
p = 0.
tspan = (0.,tf)
```
where `de_chan0` is the initial condition for the delay channel, which is a vector of arrays whose `k`th entry stores the schduled delay time for `k`th delay channel. Here we assume $X_I(0) = 0$, thus only an empty array.
Next, we choose a delay SSA algorithm `DelayDirect()` and define the problem

```julia
aggregatoralgo = DelayDirect()
save_positions = (false,false)
dprob = DiscreteProblem(u0, tspan, p)
jprob = JumpProblem(dprob, aggregatoralgo, jumpset, save_positions = (false,false))
djprob = DelayJumpProblem(jprob,delaysets,de_chan0)
```
where `DelayJumpProblem` inputs `JumpProblem`, `DelayJumpSet` and the initial condition of the delay channel `de_chan0`.

## Visualisation

Now we can solve the problem and plot a trajectory

```julia
sol = solve(djprob, SSAStepper(), seed=2, saveat =.1, save_delay_channel = false)
```

Then we simulate $10^4$ trajectories and calculate the evolution of mean value for each reactant

```julia
using DiffEqBase
Sample_size = Int(1e4)
ens_prob = EnsembleProblem(djprob)
ens =@time solve(ens_prob,SSAStepper(),EnsembleThreads(),trajectories = Sample_size, saveat = .1, save_delay_channel =false)
```
![degradation1](../assets/delay_degradation1.svg)

### Verification with the exact solution
Lastkt, we can compare with the mean values of the exact solutions $X_I, X_A$
```julia
timestamps = 0:0.1:tf
a = β + γ 
mean_x_A(t) = C/a*(1-exp(-a*t))
mean_x_I(t)= 0<=t<=τ ? C*β/(a-γ)*((1-exp(-γ*t))/γ - (1-exp(-a*t))/a) : C*β/a*((1-exp(-γ*τ))/γ + exp(-a*t)*(1-exp((a-γ)τ))/(a-γ))
```
![degradation2](../assets/delay_degradation2.svg)
![degradation3](../assets/delay_degradation3.svg)



## Reference
[1] Lafuerza, L. F., & Toral, R. (2011). Exact solution of a stochastic protein dynamics model with delayed degradation. Physical Review E, 84(5), 051121.