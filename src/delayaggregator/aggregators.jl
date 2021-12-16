abstract type AbstractDelayAggregatorAlgorithm end

struct DelayRejection <: AbstractDelayAggregatorAlgorithm end
struct DelayDirect <: AbstractDelayAggregatorAlgorithm end
struct DelayDirectCR <: AbstractDelayAggregatorAlgorithm end
struct DelayMNRM <: AbstractDelayAggregatorAlgorithm end

is_spatial(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_vartojumps_map(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::AbstractDelayAggregatorAlgorithm) = false
needs_depgraph(aggregator::DelayMNRM) = true
needs_depgraph(aggregator::DelayDirectCR) = true