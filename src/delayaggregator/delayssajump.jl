"""
An aggregator interface for Delay SSA-like algorithms.

### Required Fields
- `next_jump`          # the next jump to execute
- `prev_jump`          # the previous jump that was executed
- `next_jump_time`     # the time of the next jump
- `end_time`           # the time to stop a simulation
- `cur_rates`          # vector of current propensity values
- `sum_rate`           # sum of current propensity values
- `ma_jumps`           # any MassActionJumps for the system (scalar form)
- `rates`              # vector of rate functions for ConstantRateJumps
- `affects!`           # vector of affect functions for ConstantRateJumps
- `save_positions`     # tuple for whether to save the jumps before and/or after event
- `rng`                # random number generator
- `next_delay`         # the index of the delay channel 
- `num_next_delay`     # how many times needed for updating the states in the next_delay channel
- `time_to_next_jump`  # the time to the next jump (time gap)
- `dt_delay`           # the time to the next delay reaction

### Optional fields:
- `dep_gr`             # dependency graph, dep_gr[i] = indices of reactions that should
                       # be updated when rx i occurs.    
"""
abstract type AbstractDSSAJumpAggregator <: AbstractJumpAggregator end

DiscreteCallback(c::AbstractDSSAJumpAggregator) = DiscreteCallback(c, c, initialize = c, save_positions = c.save_positions)

########### The following routines are templates for all Delay SSAs ###########
########### Generally they should not need to be overloaded.  ###########

## Users will normally define:
# aggregate
# initialize!
# execute_jumps!
# generate_jumps!

# condition for jump to occur
@inline function (p::AbstractDSSAJumpAggregator)(u, t, integrator)
    p.next_jump_time == t
end

# executing jump at the next jump time
function (p::AbstractDSSAJumpAggregator)(integrator)
    execute_jumps!(p, integrator, integrator.u, integrator.p, integrator.t)
    generate_jumps!(p, integrator, integrator.u, integrator.p, integrator.t)
    register_next_jump_time!(integrator, p, integrator.t)
    nothing
end

# setting up a new simulation for initialize
function (p::AbstractDSSAJumpAggregator)(dj, u, t, integrator) # initialize
    initialize!(p, integrator, u, integrator.p, t)
    register_next_jump_time!(integrator, p, integrator.t)
    u_modified!(integrator,false)
    nothing
end

############################## Generic Routines ###############################


"""
    fill_rates_and_sum!(p::AbstractDSSAJumpAggregator, u, params, t)

Reevaluate all rates and their sum.
"""
function fill_rates_and_sum!(p::AbstractDSSAJumpAggregator, u, params, t)
    sum_rate = zero(typeof(p.sum_rate))

    # mass action jumps
    majumps   = p.ma_jumps
    cur_rates = p.cur_rates
    @inbounds for i in 1:get_num_majumps(majumps)
        cur_rates[i] = evalrxrate(u, i, majumps)
        sum_rate    += cur_rates[i]
    end

    # constant rates
    rates = p.rates
    idx   = get_num_majumps(majumps) + 1
    @inbounds for rate in rates
        cur_rates[idx] = rate(u, params, t)
        sum_rate += cur_rates[idx]
        idx += 1
    end

    p.sum_rate = sum_rate
    nothing
end


"""
    calculate_jump_rate(ma_jumps, rates, u, params, t, rx)

Recalculate the rate for the jump with index `rx`.
"""
@inline function calculate_jump_rate(ma_jumps, num_majumps, rates, u, params, t, rx)
    if rx <= num_majumps
        return evalrxrate(u, rx, ma_jumps)
    else
        @inbounds return rates[rx - num_majumps](u, params, t)
    end
end


@inline function dt_delay_generation!(p, integrator)
    next_jump = p.next_jump
    ttnj = p.time_to_next_jump
    @unpack delay_trigger_set, delay_interrupt_set = integrator.delayjumpsets
    if p.next_delay != nothing || any(set->next_jump in set, [delay_interrupt_set, delay_trigger_set])
        find_next_delay_dt!(p, integrator)
    else
        p.dt_delay -= ttnj
    end
end

@inline function update_state_delay!(p::AbstractDSSAJumpAggregator, integrator, u, t)
    @unpack ma_jumps, next_jump, next_delay, num_next_delay, next_jump_time, time_to_next_jump = p
    @unpack delay_interrupt, delay_trigger, delay_trigger_set, delay_interrupt_set, delay_complete = integrator.delayjumpsets

    ttnj = time_to_next_jump
    if next_delay == nothing
        # update_state ! 
        num_ma_rates = get_num_majumps(ma_jumps)
        if next_jump <= num_ma_rates # is next jump a mass action jump
            if u isa SVector
                integrator.u = executerx(u, next_jump, ma_jumps)
            else
                @inbounds executerx!(u, next_jump, ma_jumps)
            end
        else
            idx = next_jump - num_ma_rates
            @inbounds p.affects![idx](integrator)
        end
        # shift delay channel !
        shift_delay_channel!(integrator.de_chan,ttnj) #shift delays in all channel
        if next_jump in delay_interrupt_set
            update_delay_interrupt!(p, integrator)
        elseif next_jump in delay_trigger_set
            update_delay_trigger!(p, integrator)
        end
    else
        shift_delay_channel!(integrator.de_chan, ttnj)
        update_delay_channel!(integrator.de_chan)
        update_delay_complete!(p, integrator)
    end
    # save jump that was just executed
    p.prev_jump = next_jump
    nothing
end

"""
    function compare_delay!(p::AbstractDSSAJumpAggregator, de_chan, dt_delay, dt_reaction, t)

Compare dt_delay and dt_reaction.
"""
function compare_delay!(p::AbstractDSSAJumpAggregator, de_chan, dt_delay, dt_reaction, t)
    if  dt_reaction < dt_delay
        ttnj = dt_reaction
        next_delay = nothing
        num_next_delay = nothing
    else
        ttnj = dt_delay
        next_delay, num_next_delay = find_next_delay_vec(de_chan, ttnj)
    end
    p.time_to_next_jump = ttnj
    p.next_delay = next_delay
    p.num_next_delay = num_next_delay
    @fastmath p.next_jump_time = t + ttnj
    nothing
end

"""
    function find_next_delay_dt!(p, integrator)

Find the minimal dt_delay in all delay channels.
"""
function find_next_delay_dt!(p, integrator)
    de_chan = integrator.de_chan
    T = typeof(integrator.t)
    val_vec = Vector{T}(undef,length(de_chan))
    @inbounds for i in eachindex(de_chan)
        val_vec[i] = isempty(de_chan[i]) ? typemax(T) : minimum(de_chan[i])
    end
    p.dt_delay = minimum(val_vec)
    nothing
end

function find_next_delay_num!(p, de_chan::Vector{Vector{T}}) where {T}
    @label restart
    val_vec = Vector{T}(undef,length(de_chan))
    @inbounds for i in eachindex(de_chan)
        val_vec[i] = isempty(de_chan[i]) ? typemax(T) : minimum(de_chan[i])
    end
    val = minimum(val_vec)
    if val < eps(T)
        shift_delay_channel!(de_chan, eps(T))
        update_delay_channel!(de_chan)
        @goto restart
    end
    p.next_delay, p.num_next_delay = find_next_delay_vec(de_chan, val)
    p.time_to_next_jump = val
    nothing
    # return val
end

@inline function shift_delay_channel!(de_chan::Vector{Vector{T1}},ttnj::T2) where {T1<:Real,T2<:Real}
    @inbounds for idx in eachindex(de_chan)
        for j in eachindex(de_chan[idx])
            de_chan[idx][j] -= ttnj
        end
    end
end

@inline function update_delay_channel!(de_chan::Vector{Vector{T}}) where {T<:Real}
    for idx in eachindex(de_chan)
        filter!(x->x>zero(T), de_chan[idx])
    end
end
"""
    function update_delay_interrupt!(p, integrator)
"""
function update_delay_interrupt!(p, integrator)
    delay_interrupt = integrator.delayjumpsets.delay_interrupt
    execute_delay_interrupt!(delay_interrupt[p.next_jump], p, integrator)
end
function execute_delay_interrupt!(delay_interrupt_affect!::Function, p, integrator)
    delay_interrupt_affect!(integrator, p.rng)
end

"""
    function update_delay_trigger!(p, integrator)
"""
function update_delay_trigger!(p, integrator)
    delay_trigger = integrator.delayjumpsets.delay_trigger

    execute_delay_trigger!(delay_trigger[p.next_jump], p, integrator)
end
function execute_delay_trigger!(delay_trigger_affect!::Function, p, integrator)
    delay_trigger_affect!(integrator, p.rng)
end
function execute_delay_trigger!(delay_trigger_affect!::Vector{Pair{Int64,T}}, p, integrator) where {T}
    for (chan_idx, τ) in delay_trigger_affect!
        append!(integrator.de_chan[chan_idx], τ)
    end
end

"""
    function update_delay_complete!(p, integrator)

This function modifies `integrator.u` and `integrator.de_chan` upon delay completion
"""
function update_delay_complete!(p, integrator)
    delay_complete = integrator.delayjumpsets.delay_complete
    next_delay_vec, num_next_delay_vec = p.next_delay, p.num_next_delay
    @inbounds for j in eachindex(next_delay_vec)
        next_delay = next_delay_vec[j]
        num_next_delay = num_next_delay_vec[j]
        execute_delay_complete!(delay_complete[next_delay], num_next_delay, integrator, p.rng)
    end
end

function execute_delay_complete!(delay_complete::Vector{Pair{Int64,Int64}}, num_next_delay::Int64, integrator, rng)
    @inbounds for (i, ξ) in delay_complete
        if integrator.u isa SVector
            # integrator.u = setindex(integrator.u, integrator.u[i] + num_next_delay*ξ, i)
            integrator.u = setindex!!(integrator.u, integrator.u[i] + num_next_delay*ξ, i)
        else     
            integrator.u[i] += num_next_delay*ξ
        end
    end
end
function execute_delay_complete!(delay_complete::Function, num_next_delay::Int64, integrator, rng)
    for _ in 1:num_next_delay
        delay_complete(integrator, rng)
    end
end

"""
    function dep_gr_delay(p, integrator)

Generate delay dependency graph
input::Int next_delay idx
output::Dict next_delay idx => reactions need to be updated
"""
function dep_gr_delay(p::AbstractDSSAJumpAggregator, integrator)
    # num_delay_chan = length(integrator.de_chan)
    dict_ = Dict{Int,Vector{Int}}()
    dict_complete = integrator.delayjumpsets.delay_complete
    # var_to_jumps = var_to_jumps_map(length(integrator.u),p.ma_jumps)
    @inbounds for key in keys(dict_complete)
        delay_complete_action = dict_complete[key]
        if typeof(delay_complete_action)<:Vector{Pair{Int,Int}}
            vars = first.(delay_complete_action)
        else
            vars = vec(1:length(integrator.u))
        end
        jumps = unique(vars_to_jumps_delay(p, vars))
        push!(dict_, key=>jumps)
    end
    dict_
end
function vars_to_jumps_delay(p, vars)
    jumps = []
    for i in eachindex(p.dep_gr)
        for var in vars
            if var in p.dep_gr[i]
                push!(jumps, i)
            end
        end
    end
    return jumps
end

"""
    function find_num_in_vec(A::Vector{Vector{T}}, position_index::Vector{Int64}, x::T)

Find the number of values which in each vector elements equal to `x` according to the corresponding index position specified by the element in the `position_index` vector in the given vetcer `A`.

# Examples
```julia-repl
julia> A =  [[0.09,0.09,0.1],[0.3,0.09,0.1],[0.09]]
3-element Vector{Vector{Float64}}:
 [0.09, 0.09, 0.1]
 [0.3, 0.09, 0.1]
 [0.09]

julia> position_index =  [1,2,3]
3-element Vector{Int64}:
 1
 2
 3

julia> find_num_in_vec(A::Vector, position_index::Vector{Int64}, 0.09)
3-element Vector{Int64}:
2
1
1
```
"""
@inline function find_num_in_vec(A::Vector{Vector{T}}, position_index::Vector{Int64}, x::T) where {T}
    number_in_vec = Vector{Int64}(undef, length(position_index))
    @inbounds for i in eachindex(position_index)
        number_in_vec[i] = count(==(x),A[position_index[i]])
    end
    return number_in_vec
end

"""
    find_next_delay_vec(A::Vector{Vector{T}}, x::T)

Returns two vectors. The first is the position index vector `position_index` of vector `A`, and the second is the vector `num_in_vec` composed of the number of values in the position index corresponding to `position_index` equal to `x`.

find the minimal dt_delay in various delay channel.

# Examples
```julia-repl
julia> A =  [[0.09,0.09,0.1],[0.3,0.09,0.1],[0.09]]
3-element Vector{Vector{Float64}}:
 [0.09, 0.09, 0.1]
 [0.3, 0.09, 0.1]
 [0.09]

 julia> x=0.09
 0.09

julia> find_next_delay_vec(A, x)
([1, 2, 3], [2, 1, 1])
```
"""
function find_next_delay_vec(de_chan::Vector{Vector{T}}, ttnj::T) where {T}
    position_index = findall(de_chan->ttnj in de_chan, de_chan)
    num_in_vec = find_num_in_vec(de_chan, position_index, ttnj)
    return position_index, num_in_vec
end
