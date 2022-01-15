# Delay Direct Method Algorithm 

The number of discarded $\Delta$ will be approximately equal to the number of delayed reactions that initiate. This follows because, other than the stored completions at the time the code terminates, every delayed completion will cause one computed $\Delta$ to be discarded. Thus, Cai [1] developed an algorithm, called the Direct Method for systems with delays, in which no random variables are discarded.

The principle of Direct Method is the same as that of the original Gillespie Algorithm and the Rejection Method above: use one random variable to calculate when the next reaction initiates and use another random variable to calculate which reaction occurs at that future time. However, Direct Method updates the state of the system and propensity functions due to stored delayed reactions during the search for the next initiation time. In this way it is ensures that no random variables are discarded as in the Rejection Method. 

## Algorithm

Suppose that at time $t$ there are ongoing delayed reactions set to complete at times $t+T_1, t+T_2, \ldots, t+T_d$. Define $T_0=0$ and $T_{d+1}=\infty$.

Define *Tstruct*, whose $i$th $(i=1,\dots,d)$ row stores $T_i$ and the index, $\mu_i$, of the reaction that $T_i$ is associated with.
1. Initialize. Set the initial number of molecules of each species and set  $t=0$. Clear *Tstruct*.
2. Calculate the propensity of function $a_k$, for each reaction $k \in 1,\ldots, M$.
3. Set $a_0=\sum_{k=1}^M{a_k}$.
4. Generate  $\Delta$.
   - Input the time $t$ and $a_0=\sum_{k=1}^M{a_k}$.
   - Generate an independent uniform $(0,1)$ random number $r_1$.
   - If *Tstruct* is empty 
     - This means there is no ongoing delay reactions, set $\Delta = 1/a_0\ln(1/r_1)$.
   - Else
     - Set $i=1$, $a_t = a_0T_1$ and  $F=1-e^{-a_t}$.
     -  While $F < r_1$
         - Update the state vector $\mathbf{x}$ due to the finish of the delayed reaction $t+T_i$.
         - Calculate propensity $a_k(t+T_{i+1})$ due to the finish of the delayed reaction at $t+T_{i+1}$ and calculate $a_0(t+T_{i+1})$.
         - Update $a_t=a_t+a_0(t+T_{i+1})(T_{i+1}-T_i)$.
         - Update $F=1-e^{-a_t},i=i+1$.
     - EndWhile
     - Calculate Calculate propensity $a_k(t+T_i)$ due to the finish of the delayed reaction at $t+T_i$ and calculate $a_0(t+T_i)$.
     - Set $\Delta=T_i-(\ln(1-r_1)+a_t-a_0(t+T_i)(T_{i+1}-T_i))/a_0(t+T_i)$.
   - EndIf
5. If $\Delta\in[T_i,T_{i+1})$, delete the columns $1,\ldots,i$ of $T_i$ and set $T_j=T_j-\Delta$.
6. Generate an independent uniform $(0,1)$ random number $r_2$.
7. Find $\mu\in[1,\dots,m]$ such that
   ```math
   \sum_{k=1}^{\mu-1} a_k < r_2 \leq \sum_{k=1}^{\mu}a_k
   ```
   where the $a_k$ and $a_0$ are generated in step 4.
8. If $\mu\in \text{ND}$ , update the number of each molecular species according to the reaction $\mu$.
9. If $\mu\in \text{CD}$, update *Tstruct* by adding the row $[\tau_\mu,\mu]$ so that $Tstruct(i,1)<Tstruct(i+1,1)$ still holds for all **i**.
10. If $\mu\in \text{ICD}$, update the system according to the initiation of $\mu$ and update *Tstruct* by adding the row $[\tau_\mu,\mu]$ so that $Tstruct(i,1)<Tstruct(i+1,1)$ still holds for all $i$.
11. Set $t=t+\Delta$.
12. Return to Step 2 or quit.

Remark. Notice that in the above pseudo-code, we modified the Step 4 in the orignal algorithm for computational efficiency but both are equivalent.

## References

[1] Xiaodong Cai, "Exact stochastic simulation of coupled chemical reactions with delays", The Journal of Chemical Physics 126, 124108(2007).
[https://doi/10.1063/1.2710253](https://aip.scitation.org/doi/10.1063/1.2710253)