"""
$(TYPEDEF)
An abstract type that contains delay stochastic simulation algorithms:
- DelayDirect
- DelayRejection
- DelayMNRM
- DelayDirectCR
"""
abstract type AbstractDelayAggregatorAlgorithm end
"""
$(TYPEDEF)
Delay Rejection Method from Barrio, Manuel, Kevin Burrage, André Leier, and Tianhai Tian. "Oscillatory regulation of Hes1: discrete stochastic delay modelling and simulation." PLoS computational biology 2, no. 9 (2006): e117..
"""
struct DelayRejection <: AbstractDelayAggregatorAlgorithm end
"""
$(TYPEDEF)
Delay Direct Method from  Xiaodong Cai, "Exact stochastic simulation of coupled chemical reactions with delays", The Journal of Chemical Physics 126, 124108(2007).
"""
struct DelayDirect <: AbstractDelayAggregatorAlgorithm end
"""
$(TYPEDEF)
A modifed Composition-Rejection Direct Method with delays (DelayDirectCR), implementation combining features from the original article and from the code in  `JumpProcesses` package: DirectCR :
*A constant-time kinetic Monte Carlo algorithm for simulation of large biochemical reaction networks*,
by A. Slepoy, A.P. Thompson and S.J. Plimpton, J. Chem. Phys, 128, 205101 (2008).
and
*Efficient Formulations for Exact Stochastic Simulation of Chemical Systems*,
by S. Mauch and M. Stalzer, ACM Trans. Comp. Biol. and Bioinf., 8, No. 1, 27-35 (2010).
"""
struct DelayDirectCR <: AbstractDelayAggregatorAlgorithm end
"""
$(TYPEDEF)
A modifed version of the Delay Next Reaction Method from
David F. Anderson, "A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays", The Journal of Chemical Physics 128, 109903(2008).
"""
struct DelayMNRM <: AbstractDelayAggregatorAlgorithm end
"""
$(TYPEDEF)
An adaptaton of the COEVOLVE algorithm for simulating any compound jump process
that evolves through time. This method handles variable intensity rates with
user-defined bounds and inter-dependent processes. It reduces to NRM when rates
are constant.
M. Farajtabar, Y. Wang, M. Gomez-Rodriguez, S. Li, H. Zha, and L. Song,
COEVOLVE: a joint point process model for information diffusion and network
evolution, Journal of Machine Learning Research 18(1), 1305–1353 (2017). doi:
10.5555/3122009.3122050.
This is a further adaptation from https://github.com/SciML/JumpProcesses.jl/pull/276 with time delays implementation.
"""
struct DelayCoevolve <: AbstractDelayAggregatorAlgorithm end

needs_vartojumps_map(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::AbstractDelayAggregatorAlgorithm) = false

needs_depgraph(aggregator::DelayMNRM) = true
needs_vartojumps_map(aggregator::DelayMNRM) = true
needs_depgraph(aggregator::DelayDirectCR) = true
needs_vartojumps_map(aggregator::DelayDirectCR) = true
needs_depgraph(aggregator::DelayCoevolve) = true
needs_vartojumps_map(aggregator::DelayCoevolve) = true
supports_variablerates(aggregator::AbstractDelayAggregatorAlgorithm) = false
supports_variablerates(aggregator::DelayCoevolve) = true