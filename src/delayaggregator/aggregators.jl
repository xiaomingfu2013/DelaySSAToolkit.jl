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
A modifed Composition-Rejection Direct Method with delays (DelayDirectCR), implementation combining features from the original article and from the code in  `DiffEqJump` package: DirectCR :
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

needs_vartojumps_map(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::AbstractDelayAggregatorAlgorithm) = false

needs_depgraph(aggregator::DelayMNRM) = true
needs_vartojumps_map(aggregator::DelayMNRM) = true
needs_depgraph(aggregator::DelayDirectCR) = true
needs_vartojumps_map(aggregator::DelayDirectCR) = true