# Delay Rejection Method Algorithm
  Simulation methods for systems with delays need to calculate when reactions initiate and store when they complete. However, because of the delayed reactions, the propensity functions can change between initiation times. Bratsun et al. [1] and Barrio et al. [2] used an algorithm for computing the initiation times that is exactly like the original Gillespie Algorithm except that if there is a stored delayed reaction set to finish within a computed timestep, then the computed timestep is discarded, and the system is updated to incorporate the stored delayed reaction. The algorithm then attempts another step starting at its new state. This algorithm is called Rejection Method.

### Pseudo code

1. Initialize. Set the initial number of molecules of each species and set $t = 0$.
2. Calculate the propensity function, $a_k$, for each reaction.
3. Set $a_0 = \sum_{k=1}^M a_k$.
4. Generate an independent uniform $(0,1)$ random number, $r_1$, and set $\Delta = 1/a_0 \ln(1/r_1)$.
5. If there is a delayed reaction set to finish in $[t, t + \Delta)$
   - Discard $\Delta$.
   - Updatetto be the time of the next delayed reaction,$\mu$.
   - Updatexaccording to the stored reaction $\mu$.
6. Else
   - Generate an independent uniform$(0,1)$ random number $r_2$.
   - Find $\mu\in[1,\ldots, m]$ such that
   ```math
   \sum_{k=1}^{\mu-1} a_k(t) < r_2 a_0 < \sum_{k=1}^\mu a_k(t)
   ```
   - If $\mu\in \text{ND}$, update the number of each molecular species according to reaction $\mu$.
   - If $\mu\in \text{CD}$, store the information that at time $t+\tau_\mu$ the system must be updated according to reaction $\mu$.
   - If $\mu\in \text{ICD}$, update the system according to the initiation of $\mu$ and store that at time $t+\tau_\mu$ the system must be updated according to the completion of reaction $\mu$.
   - Set $t = t +\Delta$
7. Endif
8. Return to step 2 or quit.


[1] Dmitri A. Bratsun, Dmitri N. Volfson, Jeff Hasty, and Lev S. Tsimring "Non-Markovian processes in gene regulation (Keynote Address)", Proc. SPIE 5845, Noise in Complex Systems and Stochastic Dynamics III, (23 May 2005).
[https://doi.org/10.1117/12.609707](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/5845/1/Non-Markovian-processes-in-gene-regulation/10.1117/12.609707.full)

[2]  Manuel Barrio, Kevin Burrage, André Leier, Tianhai Tian. "Oscillatory Regulation of Hes1: Discrete Stochastic Delay Modelling and Simulation", PLoS Computational Biology, 10.1371(2006).
[https://doi.org/10.1371/journal.pcbi.0020117](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0020117)

