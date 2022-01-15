# Delay Modified Next Reaction Method Algorithm

## Representation using Poisson processes 

Let $v_k$, $v'_k\in Z^N_{\geq 0}$ be the vectors representing the number of each species consumed and created in the $k$th reaction, respectively. Then, if $N_k(t)$ is the number of initiations of reaction $k$ by time $t$, the state of the system at time $t$ is
```math
X(t)=X(0)+\sum_{k=1}^M{N_k(t)(v^{'}_k-v_k)}
```
However, based upon the fundamental premise of stochastic chemical kinetics, $N_k(t)$ is a counting process with intensity $a_k(X(t))$ such that $\text{P}(N_k(t+\Delta t)-N_k(t)=1|X(s),s\leqslant t)=a_k(X(t))\Delta t$ for small $\Delta t$. Therefore, based upon the [counting process interpretation](https://en.wikipedia.org/wiki/Poisson_point_process), we have
```math
N_k(t)=Y_k\Big(\int^t_0a_k(X(s))ds\Big),\tag{1} 
```
where the $Y_k$ are independent, unit rate Poisson processes. Thus, $X(t)$ can be represented as the solution to the following equation:
```math
X(t)=X(0)+\sum_{k=1}^M{Y_k\Big(\int^t_0a_k(X(s))ds\Big)(v'_k-v_k)}.
```
All of the randomness in the system is encapsulated in the $Y_k$'s and has therefore been separated from the state of the system. Thus, the system only changes when one of the $Y_k$'s changes. There are actually $M + 1$ relevant time frames in the system. The first time frame is the actual, or absolute time, $t$. However, each Poisson process $Y_k$ brings its own time frame. Thus, if we define $T_k(t)=\int^t_0a_k(X(s))ds$ for each $k$, then it is relevant for us to consider $Y_k(T_k(t))$. We will call $T_k(t)$ the "**internal time**" for reaction $k$.

We denote by $P_k$ the first firing time of $Y_k$, in the time frame of $Y_k$, which is strictly larger than  $T_k$. That is, $P_k=\min \left\{s>T_k:Y_k(s)>Y(T_k)\right\}$. The main idea of the following algorithm is that by Eq.(1) the value
```math
\Delta t_k=\frac{P_k-T_k}{a_k}
```
gives the amount of absolute time needed until the Poisson process $Y_k$ fires assuming that $a_k$ remains constant. $a_k$ does remain constant until the next reaction takes place. Therefore, a minimum of the different $\Delta t_k$ gives the time until the next reaction takes place.

Now we extend the above notions to the delay system. No matter whether a reaction is contained in $\text{ND}$, $\text{CD}$, or $\text{ICD}$, the number of initiations at absolute time $t$ will be given by
```math
\text{number of initiations of reaction } k \text{ by time } t = Y_k\Big(\int_{0}^{t} a_k(X(s)\, \mathrm{d}s\Big)
```
where the $Y_k$ are independent, unit rate Poisson processes.

Because the initiations are still given by the firing times of independent Poisson processes. Therefore, if $T_k$ is the current internal time of $Y_k$, $P_k$ the first internal time after $T_k$ at which $Y_k$ fires, and the propensity function for the $k$th reaction channel is given by $a_k$, then the time until the next initiation of reaction $k$(assuming no other reactions initiate or complete) is still given by $\Delta t_k= (P_k−T_k)/a_k$. The only change to the algorithm will be in keeping track and storing the delayed completions. To each delayed reaction channel we therefore assign a vector, $s_l$, that stores the completion times of that reaction in ascending order for $l=1,\ldots,L$, where $L$ is the total number of delay channels. Thus, the time until there is a change in the state of the system, be it an initiation or a completion, will be given by:
```math
\Delta = \min\{\min_{k\in \{1,\ldots,M\}}\{\Delta t_k\}, \min_{l\in\{1,\ldots,L\}}\{s_l(1) − t\}\},
```
where $t$ is the current time of the system. These ideas form the heart of the Next Reaction Method [1] for systems with delays.

## Algorithm

1. Initialize. Set the initial number of molecules of each species and set $t = 0$. For each $k \leq M$, set $P_k = 0$ and $T_k = 0$, and for each delayed reaction channel set $s_l = [\infty]$.
2. Calculate the propensity function, $a_k$, for each reaction.
3. Generate $M$ independent, uniform $(0,1)$ random numbers, $r_k$, and set $P_k = \ln(1/r_k)$.
4. Set $\Delta t_k = (P_k − T_k)/a_k$.
5. Set $\Delta = \min\{\min_{k\in \{1,\ldots,M\}}\{\Delta t_k\}, \min_{l\in\{1,\ldots,L\}}\{s_l(1) − t\}\}$.
6. Set $t = t + \Delta$.
7. If we chose the completion of the delayed reaction $\mu$:
   - Update the system based upon the completion of the reaction $\mu$.
   - Delete the first row of $S_\mu$.
8. Elseif reaction $\mu$ initiated and $\mu\in \text{ND}$
   - Update the system according to reaction $\mu$.
9. Elseif reaction $\mu$ initiated and $\mu\in \text{CD}$
   - Update $s_\mu$ by inserting $t + \tau_\mu$ into $s_\mu$ in the second to last position.
10. Elseif reaction $\mu$ initiated and $\mu\in \text{ICD}$
    - Update the system based upon the initiation of reaction $\mu$.
    - Update $s_\mu$ by inserting $t + \tau_\mu$ into $s_\mu$ in the second to last position.
11. For each k, set $T_k = T_k + a_k \Delta$.
12. If reaction $\mu$ initiated, let $r$ be uniform$(0,1)$ and set $P_\mu = P_\mu + \ln(1/r)$.
13. Recalculate the propensity functions, $a_k$.
14. Return to step 4 or quit.

## References

[1] David F. Anderson, "A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays", The Journal of Chemical Physics 128, 109903(2008).
[https://doi/10.1063/1.2799998](https://aip.scitation.org/doi/10.1063/1.2799998)