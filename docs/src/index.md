```@meta
CurrentModule = DelaySSAToolkit
```

# [DelaySSAToolkit](@id DelaySSAToolkit_doc)

Gillespie developed a stochastic simulation algorithm (SSA) to simulate stochastic dynamics of chemically reacting systems [1]. In SSA algorithm, it is assumed that all reactions occur instantaneously. While in many biochemical reactions, such as gene transcription and translation, it can take certain time to finish after the reactions are initiated [2]. Neglecting delays in certain cases may still produce acceptable results, but in some delay-sensitive cases, such as delay-induced oscillators, neglecting delays in simulation will lead to erroneous conclusions. To solve this problem, an exact SSA for chemical reaction systems with delays，Delay SSA [3-5] was proposed, based upon the same fundamental premise of stochastic kinetics used by Gillespie in the development of his SSA.


DelaySSAToolkit.jl is a tool developed on top of [DiffEqJump.jl](https://github.com/SciML/DiffEqJump.jl) which solves the stochastic simulation with delay and contains the following features:

## Features
- Various delay stochastic simulation algorithms are provided;
- Stochastic delay type is supported;
- Multiple delay channels and simultaneous delay reactions are supported;
- A cascade of delay reactions is supported (a delay reaction that incurs other delay reactions);
- Priority queue and dependency graph are integrated for high computational performance;
- Ecosystem with [Catalyst](https://github.com/SciML/Catalyst.jl), [DiffEqJump](https://github.com/SciML/DiffEqJump.jl), [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl) and more...

## Installation

DelaySSAToolkit can be installed through the Julia package manager:
```julia
]add https://github.com/palmtree2013/DelaySSAToolkit.jl
using DelaySSAToolkit
```
## References

[1] Daniel T. Gillespie, "Exact stochastic simulation of coupled chemical reactions", The Journal of Physical Chemistry 1977 81 (25), 2340-2361.
[https://doi.org/10.1021/j100540a008](https://doi.org/10.1021/j100540a008).

[2] Qingchao Jiang, Xiaoming Fu, Shifu Yan, Runlai Li, Wenli Du, Zhixing Cao, Feng Qian, Ramon Grima, "Neural network aided approximation and parameter inference of non-Markovian models of gene expression". Nature communications, (2021) 12(1), 1-12. [https://doi.org/10.1038/s41467-021-22919-1](https://doi.org/10.1038/s41467-021-22919-1)

[3] Barrio, Manuel, Kevin Burrage, André Leier, and Tianhai Tian. "Oscillatory regulation of Hes1: discrete stochastic delay modelling and simulation." PLoS computational biology 2, no. 9 (2006): e117. [https://doi.org/10.1371/journal.pcbi.0020117](https://doi.org/10.1371/journal.pcbi.0020117)

[4] Xiaodong Cai, "Exact stochastic simulation of coupled chemical reactions with delays", The Journal of Chemical Physics 126, 124108(2007).
[https://doi/10.1063/1.2710253](https://aip.scitation.org/doi/10.1063/1.2710253).

[5] David F. Anderson, "A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays", The Journal of Chemical Physics 128, 109903(2008).
[https://doi/10.1063/1.2799998](https://aip.scitation.org/doi/10.1063/1.2799998).

