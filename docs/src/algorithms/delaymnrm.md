# Delay Modified Next Reaction Method Algorithm

Because the initiations are still given by the firing times of independent Poisson processes. Therefore, if $T_k$ is the current internal time of $Y_k$, $P_k$ the first internal time after $T_k$ at which $Y_k$ fires, and the propensity function for the $k$th reaction channel is given by $a_k$, then the time until the next initiation of reaction $k$(assuming no other reactions initiate or complete) is still given by $\Delta t_k= (P_k−T_k)/a_k$. The only change to the algorithm will be in keeping track and storing the delayed completions. To each delayed reaction channel we therefore assign a vector, $s_k$, that stores the completion times of that reaction in ascending order. Thus, the time until there is a change in the state of the system, be it an initiation or a completion, will be given by:
```math
\Delta = \min\{\Delta t_k, s_k(1) − t\}
```
where $t$ is the current time of the system. These ideas form the heart of our Next Reaction Method [4] for systems with delays.

### Pseudo code

1. Initialize. Set the initial number of molecules of each species and set $t = 0$. For each $k \leq M$, set $P_k = 0$ and $T_k = 0$, and for each delayed reaction channel set $s_k = [\infty]$.
2. Calculate the propensity function, $a_k$, for each reaction.
3. Generate $M$ independent, uniform$(0,1)$ random numbers, $r_k$, and set $P_k = \ln(1/r_k)$.
4. Set $\Delta t_k = (P_k − T_k)/a_k$.
5. Set $\Delta = \min_k\{\Delta t_k, s_k(1) − t\}$.
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

[1] David F. Anderson, "A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays", The Journal of Chemical Physics 128, 109903(2008).
[https://doi/10.1063/1.2799998](https://aip.scitation.org/doi/10.1063/1.2799998).