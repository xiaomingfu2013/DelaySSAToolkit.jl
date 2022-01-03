```@meta
CurrentModule = DelaySSAdocs
```

# [DelaySSAToolkit](@id DelaySSAToolkit_doc)

A major assumption behind the majority of stochastic models of biochemical kinetics is the memoryless hypothesis, i.e., the stochastic dynamics of the reactants is only influenced by the current state of the system, which implies that the waiting times for reaction events obey exponential distributions. Gillespie developed a stochastic simulation algorithm (SSA) to simulate stochastic dynamics for such systems [1].  While this Markovian assumption considerably simplifies model analysis, it is dubious for modelling certain non-elementary reaction events that encapsulate multiple intermediate reaction steps [2]. To simulate such problems, several exact SSA methods for chemical reaction systems with delays (also known as delay SSA) were proposed [3-5]. 

DelaySSAToolkit.jl is a tool developed on top of [DiffEqJump.jl](https://github.com/SciML/DiffEqJump.jl) which solves the stochastic simulation with delay and contains the following features:

## Features
- Various exact delay stochastic simulation algorithms are integrated;
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
and you might need to run
```julia
using Pkg
Pkg.instantiate()
```
for the first time after installation.

## Recommendation  
To solve a `DelayJumpProblem`, here are few recommendations for good performance:

- Use Catalyst.jl to build your Markovian model (model without delays). For certain algorithms that need dependency graph, it will be auto-generated. Otherwise you must explicitly construct and pass these mappings in `JumpSet` (see [Jump Problems](https://diffeq.sciml.ai/stable/types/jump_types/#Jump-Problems) for details).

- For a small number of jumps, `DelayRejection` and `DelayDirect` will often perform better than other aggregators.

- For large numbers of jumps with sparse chain like structures and similar jump rates, for example continuous time random walks, `DelayDirectCR` and `DelayMNRM` often have the best performance.

## References

[1] Daniel T. Gillespie, "Exact stochastic simulation of coupled chemical reactions", The Journal of Physical Chemistry 1977 81 (25), 2340-2361.
[https://doi.org/10.1021/j100540a008](https://doi.org/10.1021/j100540a008)

[2] Qingchao Jiang, Xiaoming Fu, Shifu Yan, Runlai Li, Wenli Du, Zhixing Cao, Feng Qian, Ramon Grima, "Neural network aided approximation and parameter inference of non-Markovian models of gene expression". Nature communications, (2021) 12(1), 1-12. [https://doi.org/10.1038/s41467-021-22919-1](https://doi.org/10.1038/s41467-021-22919-1)

[3] Barrio, Manuel, Kevin Burrage, André Leier, and Tianhai Tian. "Oscillatory regulation of Hes1: discrete stochastic delay modelling and simulation." PLoS computational biology 2, no. 9 (2006): e117. [https://doi.org/10.1371/journal.pcbi.0020117](https://doi.org/10.1371/journal.pcbi.0020117)

[4] Xiaodong Cai, "Exact stochastic simulation of coupled chemical reactions with delays", The Journal of Chemical Physics 126, 124108(2007).
[https://doi/10.1063/1.2710253](https://aip.scitation.org/doi/10.1063/1.2710253)

[5] David F. Anderson, "A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays", The Journal of Chemical Physics 128, 109903(2008).
[https://doi/10.1063/1.2799998](https://aip.scitation.org/doi/10.1063/1.2799998)

