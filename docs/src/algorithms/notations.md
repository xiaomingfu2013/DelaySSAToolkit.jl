# Notations and Basic Concepts
Consider a system consisting of $N≥1$ chemical species,$\{X_1, . . . , X_N\}$, undergoing $M ≥ 1$ chemical reactions through reaction channels $\{R_1,...,R_M\}$, each of which is equipped with a propensity function (or intensity function in the mathematics literature),$a_k(X)$. The dynamic state of this chemical system can be described by the state vector $X(t) =[X_1(t),...,X_N(t)]^T$, where $X_n[t],n = 1,...,N,$ is the number of $X_n$ molecules at time $t$, and $[·]^T$ denotes the transpose of the vector in the bracket.

  Delays, $\tau_k > 0$, in systems are between the initiation and completion of some, or all, of the reactions. Notice that the definition of $\tau_k$  is not the next reaction time of the Next Reaction Method. We partition the reactions into three sets, those with no delays, denoted ND, those that change the state of the system only upon completion, denoted CD, and those that change the state of the system at both initiation and completion, denoted ICD. The following assumption is based upon physical principles and serves as the base assumption for simulation methods of chemically reacting systems with delays:

```math
\begin{aligned}
a_k(X(t)) \Delta t + \omicron (t) = & \text{the probability that  reaction }k \\
& \text{takes place in a small time interval }[t, t + \Delta t)
\end{aligned}
```

where $\omicron (\Delta t)/\Delta t \rightarrow 0$  as  $\Delta t \rightarrow 0$.

  Thus, no matter whether a reaction is contained in ND, CD, or ICD, the number ofinitiationsat absolute timetwill be given by

```math
\text{number of initiations of reaction } k\text{ by time } t = Y_k(\int_{0}^{t} a_k(X(s))\, \text{d}s)
```

where the $Y_k$ are independent, unit rate Poisson processes.

  Because the assumption above, and hence equation $t$, only pertains to the initiation times of reactions we must handle the completions separately. There are three different types of reactions, so there are three cases that need consideration.

**Case 1**: If reaction $k$ is in ND and initiates at time $t$, then the system is updated by losing the reactant species and gaining the product species at the time of initiation.

**Case 2**: If reaction $k$ is in CD and initiates at time $t$, then the system is updated only at the time of completion, $t + \tau_k$, by losing the reactant species and gaining the product species.

**Case 3**: If reaction $k$ is in ICD and initiates at time $t$, then the system is updated by losing the reactant species at the time of initiation, $t$, and is updated by gaining the product species at the time of completion,$t + \tau_k$.