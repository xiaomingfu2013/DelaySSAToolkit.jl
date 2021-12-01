# DelaySSAToolkit

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://palmtree2013.github.io/DelaySSAToolkit.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://palmtree2013.github.io/DelaySSAToolkit.jl/dev)
<!-- [![Build Status](https://github.com/palmtree2013/DelaySSAToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/palmtree2013/DelaySSAToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain) -->
<!-- [![Coverage](https://codecov.io/gh/palmtree2013/DelaySSAToolkit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/palmtree2013/DelaySSAToolkit.jl) -->



`DelaySSAToolkit.jl` is a tool developed on top of `DiffEqJump.jl` which solves the stochastic simulation with **delay** and contains the following features:

## Features
- Various delay stochastic simulation algorithms are provided;
- Stochastic delay type is supported;
- Multiple delay channels and simultaneous delay reactions are supported;
- Priority queue and dependency graph are integrated for high computational performance;
- Ecosystem with Catalyst, DiffEqJump, DifferentialEquations and more...

## Installation
DelaySSAToolkit can be installed through the Julia package manager:
```julia
]add https://github.com/palmtree2013/DelaySSAToolkit.jl
using DelaySSAToolkit
```
More information is available in the [documentation](https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/). Please feel free to open issues and submit pull requests!

## Examples
### SEIR model

```julia
using DelaySSAToolkit, Catalyst
using DiffEqJump

rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r

u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2]
τ = 20.
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!) # Add τ to delay channel
delay_complete = Dict(1=>[2=>1, 3=>-1]) # Transfer from E to I after the completed delay reaction
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
jprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(true,true))
sol = solve(jprob, SSAStepper())
```
![seir](docs/src/assets/seir.pdf)

### A bursty model
Check the [`bursty`](@ref) for more details.
```julia
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
jprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(false,false))
ensprob = EnsembleProblem(jprob)
@time ens = solve(ensprob, SSAStepper(), EnsembleThreads(),saveat=timestamp, trajectories=10^5)
```
![bursty](docs/src/assets/bursty.pdf)

## References
[1]: Daniel T. Gillespie "Exact stochastic simulation of coupled chemical reactions", J. Phys. Chem. 1977, 81, 25, 2340–2361.
[https://doi.org/10.1021/j100540a008](https://pubs.acs.org/doi/10.1021/j100540a008)

[2]: Xiaodong Cai, "Exact stochastic simulation of coupled chemical reactions with delays", The Journal of Chemical Physics 126, 124108(2007).
[https://doi/10.1063/1.2710253](https://aip.scitation.org/doi/10.1063/1.2710253).
