abstract type AbstractDelayAggregatorAlgorithm end
"""
Gillespie, Daniel T. (1976). A General Method for Numerically Simulating the
Stochastic Time Evolution of Coupled Chemical Reactions. Journal of
Computational Physics. 22 (4): 403–434. doi:10.1016/0021-9991(76)90041-3.
"""
struct DelayRejection <: AbstractDelayAggregatorAlgorithm end
struct DelayDirect <: AbstractDelayAggregatorAlgorithm end
struct DelayDirectCR <: AbstractDelayAggregatorAlgorithm end
struct DelayMNRM <: AbstractDelayAggregatorAlgorithm end

is_spatial(aggregator::AbstractDelayAggregatorAlgorithm) = false

needs_vartojumps_map(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::DelayMNRM) = true
needs_depgraph(aggregator::DelayDirectCR) = true