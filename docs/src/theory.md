# Theory
## Exact Stochastic Simulation Algorithm (SSA) Without Delays

Consider a system consisting of $N≥1$ chemical species, $\{X_1, . . . , X_N\}$, undergoing $M ≥ 1$ chemical reactions through reaction channels $\{R_1,...,R_M\}$, each of which is equipped with a propensity function (or intensity function in the mathematics literature), $a_k(X)$. The dynamic state of this chemical system can be described by the state vector $X(t) =[X_1(t),...,X_N(t)]^T$, where $X_n[t],n = 1,...,N,$ is the number of $X_n$ molecules at time $t$, and $[·]^T$ denotes the transpose of the vector in the bracket.

  Following Gillespie, A the dynamics of reaction $R_m$ defined by a state-change vector $\nu_m = [\nu_{1m} ,...,\nu_{Nm}]^T$, where $\nu_{nm}$ gives the changes in the $X_n$ molecular population produced by one $R_m$ reaction, and a propensity function $a_m(t)$ together with the fundamental premise of stochastic chemical kinetics:

```math
\begin{equation}
\begin{aligned}
a_m(t)dt =& \text{the probability, given } X(t)=x, \\
& \text{that one reaction }R_m \text{ will occur in the}\\
& \text{next infinitesimal time interval }[t,t+d_t].
\end{aligned}
\end{equation}
```

  Defining the probability rate constant $c_m$ as the probability that a randomly selected combination of $R_m$ reactant molecules reacts in a unit time period, we can calculate  $a_m(t)$ fromcmand the molecular numbers ofRmreactants at time $t$ using the method given by Gillespie.

  For a chemical system in a given state $X(t)=x$ at time $t$,assuming that all reactions occur instantly, Gillespie’s exact SSA answers the following two questions: (i)  when will the next reaction occur?  (ii)  which reaction will occur? Specifically, Gillespie’s exact SSA simulates the following event in each step:

```math
\begin{equation}
\begin{aligned}
\text{E:} & \text{no reaction occurs in the time interval }[t,t+\tau],\\
& \text{and a reaction }R_\mu \ \text{occurs in the infinitesimal}\\
& \text{time interval }[t+\tau,t+\tau+d_\tau].
\end{aligned}
\end{equation}
```

  Based upon the fundamental premise Eq. (1), Gillespie showed that that $\tau$ and $\mu$ are two independent random variables and have the following probability density functions (PDFs), respectively:

```math
\begin{equation}
f_\tau(\tau)=a_0(t) \exp(-a_0(t)\tau), \tau > 0,
\end{equation}
```

and

```math
\begin{equation}
f_\mu(\mu)=a_\mu(t)/a_0(t), \mu = 1,...,M,
\end{equation}
```

where $a_0(t)=\begin{matrix} \sum_{m=1}^M a_m(t) \end{matrix}$. According to the PDF Eq.(4), a realization of $\mu$ can be generated from a standard uniform random variable $u_1$, by taking $\mu$ to be the integer for which $\begin{matrix} \sum_{j=1}^{\mu-1} a_j(t) \end{matrix} < u_1 a_0(t) ≤ \begin{matrix} \sum_{j=1}^\mu a_j(t) \end{matrix}$;based on the PDF Eq.(3), a realization of $\tau$ can be generated from another standard uniform random variable $u_2$ as $\tau=−\ln(u_2)/a_0(t)$. Therefore, Gillespie’s exact SSA generates a realization of $\mu$ and $\tau$ in each step of simulation, and then updates the time and system state as $t\leftarrow t+\tau$ and  $\mathbf{x} \leftarrow \mathbf{x}+ \mathbf{\nu_\mu}$, respectively.

## Exact SSA For Coupled Chemical Reaction With Delays

### Delay Direct method

  As in the derivation of Gillespie’s exact SSA, we first need to find the probability of event Eq.(2), that is defined as $P(\tau,\mu)d\tau$, where $P(\tau,\mu)$ is the joint PDF of $\tau$ and $\mu$. Suppose that there are $d$ ongoing reactions at timet, which will finish at $t+T_1,...,t+T_{d}$, respectively. Without loss of generality, we assume that $T_1≤T_2≤...≤T_d$. Unlike in the reaction system without delays where the propensity functions remain unchanged in the time interval $[t,t+\tau]$, the propensity functions here change at $t+T_i,i=1,...,d$, due to delayed reactions. We need to take into account such changes in propensity functions when deriving  $P(\tau,\mu)$.

  As in the derivation of Gillespie’s exact SSA, $P(\tau,\mu)d\tau$ can be found from the fundamental premise Eq.(1) as

```math
\begin{equation}
P(\tau,\mu)d\tau=P_0(\tau) a_\mu(\tau,\mu)d\tau,
\end{equation}
```

where $P_0(\tau)$ is the probability that no reaction will occur in the time interval $[t,t+\tau]$, while $a_\mu(t+\tau)d\tau$ is the probability that a reaction $R_\mu$ occurs in $[t+\tau,t+\tau+d\tau]$. Defining $T_0=0$ and $T_{d+1}=\infty$, we can find $P_0(\tau)$ for $\tau$ that lies in different time intervals $[T_i,T_{i+1}),i=0,...,d$. If $\tau \in [T_i,T_i+1)$, we define the event $E_j$ as the event that no reaction occurs in the time interval $[t+T_j,t+T_j+1),j=0,...,j=i−1$, respectively,and the event  $E_i$  as the event that no reaction occurs in the time interval $[t+T_i,t+\tau)$. Then, we can express $P_0(\tau)$ as

```math
\begin{equation}
P_0(\tau)=P(E_0,...,E_i)=P(E_0) \prod_{j=1}^i P(E_j丨E_0,...,E_{j-1}).
\end{equation}
```

  From the derivation of Gillespie’s exact SSA,we know that  $P(E_0) = \exp (−a_0(t)T_1)$,  $P(E_j丨E_0,...,E_{j-1}) = \exp(-a_0(t+T_j)T_1) × (T_{j+1}−T_j),j=0,...,i−1$,   and   $P(E_i丨E_0,...,E_{i-1}) = \exp(-a_0(t+T_i)(\tau-T_i))$.  Notice that propensity functions change at $t+T_j$ after a delayed reaction finishes, and we use $a_0(t+T_j)$ to represent the new $a_0$. The probability $P_0(\tau)$ is then given by

```math
\begin{equation}
\begin{aligned}
& P_0(\tau) = \exp \bigg (-\sum_{j=0}^{i-1} a_0(t+T_j)(T_{j+1}-T_j)-a_0(t+T_i)(\tau-T_i) \bigg ), \\
& \tau \in [T_i,T_i+1), i = 0,...,d,
\end{aligned}
\end{equation}
```

where we assume that the first term of the exponent is equal to zero when $i = 0$. Since $P_0(\tau)$ does not depend on individual propensity functions, as shown in Eq.(7), it is seen from Eq.(5) that $\tau$ and $\mu$ are independent random variables. Combining Eq.(5) and Eq.(7) and noticing that $a_\mu(t+\tau)=a_\mu(t+T_i)$ for $\tau \in [T_i,T_i+1)$, we obtain the PDF of $\tau$ and $\mu$ as follows:

```math
\begin{equation}
\begin{aligned}
& f_\tau(\tau) = a_0(t+T_i) \exp \bigg (-\begin{matrix} \sum_{j=0}^{i-1} a_0(t+T_j)(T_{j+1}-T_j) \end{matrix} - a_0(t+T_i)(\tau-T_i) \bigg ), \\
& \tau \in [T_i,T_i+1), i = 0,...,d,
\end{aligned}
\end{equation}
```

and

```math
f_\mu(\mu)=a_\mu(t+T_i)/a_0(t+T_i), \mu = 1,...,M,\tau \in [T_i,T_i+1),
```

It is not difficult to verify that $\int_{0}^{\infty} f_\tau(\tau)\, d\tau = 1$. In simulation, $\mu$ can be generated, from a standard uniform random variable $u_1$, by taking $\mu$ to be the integer for which $\begin{matrix} \sum_{j=1}^{\mu-1} a_j(t+T_i) \end{matrix} < u_1 a_0(t+T_i) ≤ \begin{matrix} \sum_{j=1}^\mu a_j(t+T_i) \end{matrix}$, after $\tau$ is generated to be in the time interval $[T_i,T_{i+1})$.We next derive the method of generating  $\tau$ according to its PDF in Eq.(8).

  The cumulative distribution function of $\tau$can be found from Eq.(8) as

```math
\begin{aligned}
& F_\tau(\tau)=1 - \exp  \bigg (-\begin{matrix} \sum_{j=0}^{i-1} a_0 \end{matrix}(t+T_j)(T_{j+1}-T_j)-a_0(t+T_i)(\tau-T_i) \bigg ), \\
& \tau \in [T_i,T_i+1), i = 0,...,d,
\end{aligned}
```

Then, we can generate $\tau$ from a standard uniform random variable $u_2$, by taking $\tau=F_\tau^{−1}(u_2)$, where $F_\tau^{−1}(\cdot)$ represents the inverse of $F_\tau(\tau)$. More specifically, we can obtain $\tau$ as follows:

  Find $T_i$ such that  $F_\tau(T_i) ≤ u_2 ≤ F_\tau(T_{i+1})$, then calculate  $\tau$ from

```math
\begin{aligned}
& \tau = T_i + {{-\ln (1-u_2)-\begin{matrix} \sum_{j=0}^{i-1} a_0(t+T_j)(T_{j+1}-T_j) \end{matrix} } \over {a_0(t+T_j)}} \\
& \tau \in [T_i,T_i+1).
\end{aligned}
```

  Since we need $T_1,...,T_d$ to generate $\tau$ and $\mu$, we define an array of data structures, named $Tstruct$, whose $i$th $(i=1,...,d)$ cell stores $T_i$ and the index, $\mu_i$, of the reaction that $T_i$ is associated with. The reaction index $\mu_i$ is needed during the generation of $\tau$, when we update the propensity functions affected by the reaction that is delayed but finishes at $t+T_i$. During simulation, we need to generate $\tau$ and $\mu$, maintain *Tstruct*, and then update the state vector $X(t)$.

### Delay Rejection Method