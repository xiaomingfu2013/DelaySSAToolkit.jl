```@meta
CurrentModule = DelaySSAToolkit
```

# DelaySSAToolkit

Gillespie developed a stochastic simulation algorithm (SSA) to simulate stochastic dynamics of chemically reacting systems [1]. In SSA algorithm, it is assumed that all reactions occur instantaneously. While in many biochemical reactions, such as gene transcription and translation, it can take certain time to finish after the reactions are initiated [??]. Neglecting delays in certain cases may still produce acceptable results, but in some delay-sensitive cases, such as delay-induced oscillators, neglecting delays in simulation will lead to erroneous conclusions. To solve this problem, an exact SSA for chemical reaction systems with delays，Delay SSA [2, 3] was proposed, based upon the same fundamental premise of stochastic kinetics used by Gillespie in the development of his SSA.

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
## References

[1]: Daniel T. Gillespie "Exact stochastic simulation of coupled chemical reactions", J. Phys. Chem. 1977, 81, 25, 2340–2361.
[https://doi.org/10.1021/j100540a008](https://pubs.acs.org/doi/10.1021/j100540a008)

[2]: Xiaodong Cai, "Exact stochastic simulation of coupled chemical reactions with delays", The Journal of Chemical Physics 126, 124108(2007).
[https://doi/10.1063/1.2710253](https://aip.scitation.org/doi/10.1063/1.2710253).


