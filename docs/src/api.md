# Main API

```@meta
CurrentModule = DelaySSAdocs
```

```@docs
DelayJumpSet
DelayJumpProblem
```

## More about defining a `DelayJumpSet`
Further notes on `delay_trigger`, `delay_interrupt`, `delay_complete`.
### `delay_trigger` 
   delay_trigger defines a `Dict` type:

- Keys: Indices of reactions defined in the Markovian part that can trigger the delay reaction;
- Values: An update function or a `Pair` type that determines how to update the delay channel.
 

### `delay_interrupt`

### `delay_complete` 

## Types and Algorithms
```@docs
AbstractDelayAggregatorAlgorithm
DelayDirect
DelayRejection
DelayMNRM
DelayDirectCR
```