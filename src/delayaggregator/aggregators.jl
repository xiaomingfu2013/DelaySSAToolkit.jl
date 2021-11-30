"""
Gillespie, Daniel T. (1976). A General Method for Numerically Simulating the
Stochastic Time Evolution of Coupled Chemical Reactions. Journal of
Computational Physics. 22 (4): 403–434. doi:10.1016/0021-9991(76)90041-3.
"""
struct DelayRejection <: AbstractAggregatorAlgorithm end
struct DelayDirect <: AbstractAggregatorAlgorithm end
struct DelayMNRM <: AbstractAggregatorAlgorithm end

is_spatial(aggregator::AbstractAggregatorAlgorithm) = false

needs_vartojumps_map(aggregator::AbstractAggregatorAlgorithm) = false
needs_depgraph(aggregator::AbstractAggregatorAlgorithm) = false
needs_depgraph(aggregator::DelayMNRM) = true