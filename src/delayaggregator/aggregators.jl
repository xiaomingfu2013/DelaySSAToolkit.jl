"""
$(TYPEDEF)
An abstract type that contains delay stochastic simulation algorithms:
- DelayDirect
- 
"""
abstract type AbstractDelayAggregatorAlgorithm  end
"""
$(TYPEDEF)

"""
struct DelayRejection <: AbstractDelayAggregatorAlgorithm end
"""
$(TYPEDEF)

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

"""
struct DelayMNRM <: AbstractDelayAggregatorAlgorithm end

is_spatial(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_vartojumps_map(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::DelayMNRM) = true
needs_depgraph(aggregator::DelayDirectCR) = true