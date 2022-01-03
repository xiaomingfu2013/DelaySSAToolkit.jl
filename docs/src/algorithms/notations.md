# Notations and Basic Concepts

DelaySSAToolkit now supports four exact delay stochastic simulation algorithms, namely, delay rejection method `DelayRejction` [1,2], delay direct method `DelayDirect` [3], delay modified next reaction method `DelayMNRM` [4] and delay direct method with composition and rejection `DelayDirectCR` [5,6], which intend to solve reaction systems with different scales with best performance. We refer to [recommendation](../index.md) for details. We point out that delay direct method with composition and rejection achieves better computational efficiency by a more complex underlying data structure and a partition of propensity functions (see [5-6] for more details), and the fundamental algorithm structure remains the same as delay rejection method. Here we briefly present the algorithm stuctures of delay rejection method, delay direct method and delay modified next reaction method. First of all, we introduce some basic concepts for delay stochastic simulation algorithms.

Consider a system consisting of $N \geq 1$ chemical species, $\{X_1,\ldots, X_N\}$, undergoing $M \geq 1$ chemical reactions through reaction channels $\{R_1,\ldots,R_M\}$, each of which is equipped with a propensity function (or intensity function in the mathematics literature), $a_k(X)$. The dynamic state of this chemical system can be described by the state vector $X(t) =[X_1(t),\ldots,X_N(t)]^T$, where $X_n(t),n = 1,\ldots,N$, is the number of $X_n$ molecules at time $t$, and $[·]^T$ denotes the transpose of the vector in the bracket.

The delays, $\tau_k > 0$ for $k = 1,\ldots,d$, in systems are between the initiation and completion of some, or all, of the reactions. And $\tau_k$ is used to represent the delay time of the $k$th reaction in all delayed reactions. Notice that the definition of $\tau_k$ is not the next reaction time $\Delta$. We partition the reactions into three sets, those with no delays, denoted $\text{ND}$, those that change the state of the system only upon completion, denoted $\text{CD}$, and those that change the state of the system at both initiation and completion, denoted $\text{ICD}$. The following assumption, sometimes called the fundamental premise of chemical kinetics, is based upon physical principles and serves as the base assumption for simulation methods of chemically reacting systems with delays:

```math
\begin{aligned}
a_k(X(t)) \Delta t + \omicron (t) = & \text{ the probability that  reaction }k \\
& \text{ takes place in a small time interval }[t, t + \Delta t)
\end{aligned}
```

where $\omicron (\Delta t)/\Delta t \rightarrow 0$  as  $\Delta t \rightarrow 0$.

Because the assumption above only pertains to the initiation times of reactions we must handle the completions separately. There are three different types of reactions, so there are three cases that need consideration.

**Case 1**: If reaction $k$ is in $\text{ND}$ and initiates at time $t$, then the system is updated by losing the reactant species and gaining the product species at the time of initiation.

**Case 2**: If reaction $k$ is in $\text{CD}$ and initiates at time $t$, then the system is updated only at the time of completion, $t + \tau_k$, by losing the reactant species and gaining the product species.

**Case 3**: If reaction $k$ is in $\text{ICD}$ and initiates at time $t$, then the system is updated by losing the reactant species at the time of initiation, $t$, and is updated by gaining the product species at the time of completion,$t + \tau_k$.



## References
[1] Dmitri A. Bratsun, Dmitri N. Volfson, Jeff Hasty, and Lev S. Tsimring "Non-Markovian processes in gene regulation (Keynote Address)", Proc. SPIE 5845, Noise in Complex Systems and Stochastic Dynamics III, (23 May 2005).
[https://doi.org/10.1117/12.609707](https://doi.org/10.1117/12.609707)

[2] Manuel Barrio, Kevin Burrage, André Leier, Tianhai Tian. "Oscillatory Regulation of Hes1: Discrete Stochastic Delay Modelling and Simulation", PLoS Computational Biology, 10.1371(2006).
[https://doi.org/10.1371/journal.pcbi.0020117](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0020117)

[3] Xiaodong Cai, "Exact stochastic simulation of coupled chemical reactions with delays", The Journal of Chemical Physics 126, 124108(2007).
[https://doi/10.1063/1.2710253](https://aip.scitation.org/doi/10.1063/1.2710253)

[4] David F. Anderson, "A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays", The Journal of Chemical Physics 128, 109903(2008).
[https://doi/10.1063/1.2799998](https://aip.scitation.org/doi/10.1063/1.2799998)

[5] Slepoy, Alexander, Aidan P. Thompson, and Steven J. Plimpton. "A constant-time kinetic Monte Carlo algorithm for simulation of large biochemical reaction networks." The journal of chemical physics 128, no. 20 (2008): 05B618. [https://doi.org/10.1063/1.2919546](https://doi.org/10.1063/1.2919546)

[6] Mauch, Sean, and Mark Stalzer. "Efficient formulations for exact stochastic simulation of chemical systems." IEEE/ACM Transactions on Computational Biology and Bioinformatics 8, no. 1 (2009): 27-35. [https://doi.org/10.1109/TCBB.2009.47](https://doi.org/10.1109/TCBB.2009.47)