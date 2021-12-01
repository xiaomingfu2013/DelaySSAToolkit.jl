# Theory
## Exact Stochastic Simulation Algorithm (SSA) Without Delays

Consider a system consisting of $N \geq 1$ chemical species, $\{X_1,\ldots, X_N\}$, undergoing $M \geq 1$ chemical reactions through reaction channels $\{R_1,\ldots,R_M\}$, each of which is equipped with a propensity function (or intensity function in the mathematics literature), $a_k(X)$. The dynamic state of this chemical system can be described by the state vector $X(t) =[X_1(t),\ldots,X_N(t)]^T$, where $X_n(t),n = 1,\ldots,N,$ is the number of $X_n$ molecules at time $t$, and $[·]^T$ denotes the transpose of the vector in the bracket.

Following Gillespie, the dynamics of reaction $R_k$ defined by a state-change vector $\nu_k = [\nu_{1k} ,\ldots,\nu_{Nk}]^T$, where $\nu_{nk}$ gives the changes in the $X_n$ molecular population produced by one $R_k$ reaction, and a propensity function $a_k(t)$ together with the fundamental premise of stochastic chemical kinetics:

```math
\begin{equation}
\begin{aligned}
a_k(t)\Delta t = &\text{ the probability, given } X(t)=x, \\
&\quad \text{that one reaction }R_k \text{ will occur in the}\\
&\quad \text{next infinitesimal time interval }[t,t+\Delta t].
\end{aligned}
\end{equation}
```

Defining the probability rate constant $c_k$ as the probability that a randomly selected combination of $R_k$ reactant molecules reacts in a unit time period, we can calculate  $a_k(t)$ from $c_k$ and the molecular numbers of $R_k$ reactants at time $t$ using the method given by Gillespie.

For a chemical system in a given state $X(t)=x$ at time $t$,assuming that all reactions occur instantly, Gillespie’s exact SSA answers the following two questions: (i)  when will the next reaction occur?  (ii)  which reaction will occur? Specifically, Gillespie’s exact SSA simulates the following event in each step:

```math
\begin{equation}
\begin{aligned}
\text{E: } & \text{no reaction occurs in the time interval }[t,t+\Delta],\\
& \text{and a reaction }R_\mu \ \text{occurs in the infinitesimal}\\
& \text{time interval }[t+\Delta,t+\Delta+\Delta t].
\end{aligned}
\end{equation}
```

  Based upon the fundamental premise Eq. (1), Gillespie showed that that $\Delta$ and $\mu$ are two independent random variables and have the following probability density functions (PDFs), respectively:

```math
\begin{equation}
f_\Delta(\Delta)=a_0(t) \exp(-a_0(t)\Delta), \Delta > 0,
\end{equation}
```

and

```math
\begin{equation}
f_\mu(\mu)={{a_\mu(t)} \over {a_0(t)}}, \mu = 1,\ldots,M,
\end{equation}
```

where $a_0(t)=\begin{matrix} \sum_{k=1}^M a_k(t) \end{matrix}$. According to the PDF Eq.(4), a realization of $\mu$ can be generated from a standard uniform random variable $u_1$, by taking $\mu$ to be the integer for which $\sum_{k=1}^{\mu-1} a_k(t)  < u_1 a_0(t) \leq \sum_{k=1}^\mu a_k(t)$;based on the PDF Eq.(3), a realization of $\Delta$ can be generated from another standard uniform random variable $u_2$ as $\Delta=−\ln(u_2)/a_0(t)$. Therefore, Gillespie’s exact SSA generates a realization of $\mu$ and $\Delta$ in each step of simulation, and then updates the time and system state as $t\leftarrow t+\Delta$ and  $\mathbf{x} \leftarrow \mathbf{x}+ \mathbf{\nu_\mu}$, respectively.

## Exact SSA For Coupled Chemical Reaction With Delays

### Delay Direct method

  As in the derivation of Gillespie’s exact SSA, we first need to find the probability of event Eq.(2), that is defined as $P(\Delta,\mu)d\Delta$, where $P(\Delta,\mu)$ is the joint PDF of $\Delta$ and $\mu$. Suppose that there are $d$ ongoing reactions at timet, which will finish at $t+T_1,\ldots,t+T_{d}$, respectively. Without loss of generality, we assume that $T_1 \leq T_2 \leq \ldots \leq T_d$. Unlike in the reaction system without delays where the propensity functions remain unchanged in the time interval $[t,t+\Delta]$, the propensity functions here change at $t+T_i,i=1,\ldots,d$, due to delayed reactions. We need to take into account such changes in propensity functions when deriving  $P(\Delta,\mu)$.

  As in the derivation of Gillespie’s exact SSA, $P(\Delta,\mu)d\Delta$ can be found from the fundamental premise Eq.(1) as

```math
\begin{equation}
P(\Delta,\mu)d\Delta=P_0(\Delta) a_\mu(\Delta,\mu)d\Delta,
\end{equation}
```

where $P_0(\Delta)$ is the probability that no reaction will occur in the time interval $[t,t+\Delta]$, while $a_\mu(t+\Delta)d\Delta$ is the probability that a reaction $R_\mu$ occurs in $[t+\Delta,t+\Delta+d\Delta]$. Defining $T_0=0$ and $T_{d+1}=\infty$, we can find $P_0(\Delta)$ for $\Delta$ that lies in different time intervals $[T_i,T_{i+1}),i=0,\ldots,d$. If $\Delta \in [T_i,T_i+1)$, we define the event $E_k$ as the event that no reaction occurs in the time interval $[t+T_k,t+T_{k+1}),k=0,\ldots,k=i−1$, respectively,and the event  $E_i$  as the event that no reaction occurs in the time interval $[t+T_i,t+\Delta)$. Then, we can express $P_0(\Delta)$ as

```math
\begin{equation}
P_0(\Delta)=P(E_0,\ldots,E_i)=P(E_0) \prod_{k=1}^i P(E_k丨E_0,\ldots,E_{k-1}).
\end{equation}
```

From the derivation of Gillespie’s exact SSA,we know that  
```math
P(E_0) = \exp (−a_0(t)T_1)\\
P_0(E_k丨E_0,\ldots,E_{k-1}) = \exp(-a_0(t+T_k)T_1) × (T_{k+1}−T_k),k=0,\ldots,i−1,\\
P(E_i丨E_0,\ldots,E_{i-1}) = \exp(-a_0(t+T_i)(\Delta-T_i)).
```
Notice that propensity functions change at $t+T_k$ after a delayed reaction finishes, and we use $a_0(t+T_k)$ to represent the new $a_0$. The probability $P_0(\Delta)$ is then given by

```math
\begin{equation}
\begin{aligned}
& P_0(\Delta) = \exp \bigg (-\sum_{k=0}^{i-1} a_0(t+T_k)(T_{k+1}-T_k)-a_0(t+T_i)(\Delta-T_i) \bigg ), \\
& \Delta \in [T_i,T_i+1), i = 0,\ldots,d,
\end{aligned}
\end{equation}
```

where we assume that the first term of the exponent is equal to zero when $i = 0$. Since $P_0(\Delta)$ does not depend on individual propensity functions, as shown in Eq.(7), it is seen from Eq.(5) that $\Delta$ and $\mu$ are independent random variables. Combining Eq.(5) and Eq.(7) and noticing that $a_\mu(t+\Delta)=a_\mu(t+T_i)$ for $\Delta \in [T_i,T_i+1)$, we obtain the PDF of $\Delta$ and $\mu$ as follows:

```math
\begin{equation}
\begin{aligned}
& f_\Delta(\Delta) = a_0(t+T_i) \exp \bigg (-\begin{matrix} \sum_{k=0}^{i-1} a_0(t+T_k)(T_{k+1}-T_k) \end{matrix} - a_0(t+T_i)(\Delta-T_i) \bigg ), \\
& \Delta \in [T_i,T_i+1), i = 0,\ldots,d,
\end{aligned}
\end{equation}
```

and

```math
f_\mu(\mu)={ {a_\mu(t+T_i)} \over {a_0(t+T_i)} }, \mu = 1,\ldots,M,\Delta \in [T_i,T_i+1),
```

It is not difficult to verify that $\int_{0}^{\infty} f_\Delta(\Delta)\, d\Delta = 1$. In simulation, $\mu$ can be generated, from a standard uniform random variable $u_1$, by taking $\mu$ to be the integer for which $\begin{matrix} \sum_{k=1}^{\mu-1} a_k(t+T_i) \end{matrix} < u_1 a_0(t+T_i) ≤ \begin{matrix} \sum_{k=1}^\mu a_k(t+T_i) \end{matrix}$, after $\Delta$ is generated to be in the time interval $[T_i,T_{i+1})$.We next derive the method of generating  $\Delta$ according to its PDF in Eq.(8).

  The cumulative distribution function of $\Delta$can be found from Eq.(8) as

```math
\begin{aligned}
& F_\Delta(\Delta)=1 - \exp  \bigg (-\begin{matrix} \sum_{k=0}^{i-1} a_0 \end{matrix}(t+T_k)(T_{k+1}-T_k)-a_0(t+T_i)(\Delta-T_i) \bigg ), \\
& \Delta \in [T_i,T_i+1), i = 0,\ldots,d,
\end{aligned}
```

Then, we can generate $\Delta$ from a standard uniform random variable $u_2$, by taking $\Delta=F_\Delta^{−1}(u_2)$, where $F_\Delta^{−1}(\cdot)$ represents the inverse of $F_\Delta(\Delta)$. More specifically, we can obtain $\Delta$ as follows:

  Find $T_i$ such that  $F_\Delta(T_i) ≤ u_2 ≤ F_\Delta(T_{i+1})$, then calculate  $\Delta$ from

```math
\begin{aligned}
& \Delta = T_i + {{-\ln (1-u_2)-\begin{matrix} \sum_{k=0}^{i-1} a_0(t+T_k)(T_{k+1}-T_k) \end{matrix} } \over {a_0(t+T_i)}} \\
& \Delta \in [T_i,T_i+1).
\end{aligned}
```

  Since we need $T_1,\ldots,T_d$ to generate $\Delta$ and $\mu$, we define an array of data structures, named *Tstruct*, whose $i$th $(i=1,\ldots,d)$ cell stores $T_i$ and the index, $\mu_i$, of the reaction that $T_i$ is associated with. The reaction index $\mu_i$ is needed during the generation of $\Delta$, when we update the propensity functions affected by the reaction that is delayed but finishes at $t+T_i$. During simulation, we need to generate $\Delta$ and $\mu$, maintain *Tstruct*, and then update the state vector $X(t)$.