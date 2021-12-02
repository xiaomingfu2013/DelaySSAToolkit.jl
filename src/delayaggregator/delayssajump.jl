"""
An aggregator interface for SSA-like algorithms.

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

### Optional fields:
- `dep_gr`             # dependency graph, dep_gr[i] = indices of reactions that should
                       # be updated when rx i occurs.    
"""
abstract type AbstractDSSAJumpAggregator <: AbstractJumpAggregator end

DiscreteCallback(c::AbstractDSSAJumpAggregator) = DiscreteCallback(c, c, initialize = c, save_positions = c.save_positions)

########### The following routines are templates for all SSAs ###########
########### Generally they should not need to be overloaded.  ###########

## Users will normally define (see direct.jl for examples):
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

# setting up a new simulation
function (p::AbstractDSSAJumpAggregator)(dj, u, t, integrator) # initialize
    initialize!(p, integrator, u, integrator.p, t)
    register_next_jump_time!(integrator, p, integrator.t)
    u_modified!(integrator,false)
    nothing
end

############################## Generic Routines ###############################

"""
    register_next_jump_time!(integrator, p::AbstractDSSAJumpAggregator, t)

Adds a `tstop` to the integrator at the next jump time.
"""
# @inline function register_next_jump_time!(integrator, p::AbstractDSSAJumpAggregator, t)
#     if p.next_jump_time < p.end_time
#         add_tstop!(integrator, p.next_jump_time)
#     end
#     nothing
# end


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
    if next_jump !=nothing || next_jump in delay_trigger_set || next_jump in delay_interrupt_set
        p.dt_delay = find_next_delay_dt(integrator.de_chan)
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
        shift_delay_channel!(integrator.de_chan,ttnj) #更新 delay channel 里面的时间 for all channels
        if next_jump in delay_interrupt_set
            # delay_chan is changed according to affect_chan!
            delay_interrupt[next_jump](integrator.de_chan, p.rng) #affect_chan! decide how to modify the molecules in the delay channels
        elseif next_jump in delay_trigger_set
            # delay_trigger[next_jump] will affect de_channel
            delay_trigger[next_jump](integrator.de_chan, p.rng) # change to affect!
        end
    else
        # p.num_next_delay = num_next_delay
        shift_delay_channel!(integrator.de_chan,ttnj)
        update_delay_channel!(integrator.de_chan)
        update_delay_complete!(p, integrator, next_delay, num_next_delay, delay_complete)
    end
    # save jump that was just executed
    p.prev_jump = next_jump
    nothing
end
"""
Compare delay dt with reaction dt 
"""
function compare_delay!(p::AbstractDSSAJumpAggregator, de_chan, dt_delay, dt_reaction, t)
    if  dt_reaction < dt_delay
        ttnj = dt_reaction
        next_delay = nothing
        num_next_delay = nothing
    elseif dt_reaction >= dt_delay && dt_delay < Inf
        ttnj = dt_delay
        # next_delay = find_next_delay(de_chan)
        # num_next_delay = check_num_next_delay(de_chan[next_delay],dt_delay)
        next_delay, num_next_delay = find_next_delay_vec(de_chan, ttnj)
    else
        error("Infinite waiting time for next jump")
    end
    p.time_to_next_jump = ttnj
    p.next_delay = next_delay
    p.num_next_delay = num_next_delay
    @fastmath p.next_jump_time = t + ttnj
    nothing
end

"""
find the minimal dt_delay in various delay channel
"""

function find_next_delay_dt(de_chan::Vector{Vector{T}}) where {T}
    val_vec = Vector{T}(undef,length(de_chan))
    @inbounds for i in eachindex(de_chan)
        val_vec[i] = isempty(de_chan[i]) ? Inf : minimum(de_chan[i])
    end
    minimum(val_vec)
end

"""
Shift delay channel according to ttnj
"""

@inline function shift_delay_channel!(de_chan::Vector,ttnj)
    for idx in eachindex(de_chan)
        de_chan[idx] .-=ttnj
    end
end

"""
Update the delay channel 
"""
@inline function update_delay_channel!(de_chan::Vector)
    for idx in eachindex(de_chan)
        filter!(x->x.>0, de_chan[idx])
    end 
end

"""
Update the state upon delay completion
"""
function update_delay_complete_u!(u, next_delay_vec, num_next_delay_vec, delay_complete::Vector{Pair})
    @inbounds for j in eachindex(next_delay_vec)
        next_delay = next_delay_vec[j]
        @inbounds for (i, ξ) in delay_complete[next_delay]
            u[i] += num_next_delay_vec[j]*ξ
        end
    end
end
function update_delay_complete!(p, integrator, next_delay_vec, num_next_delay_vec, delay_complete::Dict{Int64,Any})
    @inbounds for j in eachindex(next_delay_vec)
        next_delay = next_delay_vec[j]
        num_next_delay = num_next_delay_vec[j]
        execute_delay_complete!(delay_complete[next_delay], num_next_delay, integrator, p.rng)
    end
end

function execute_delay_complete!(delay_complete::Vector{Pair{Int64,Int64}}, num_next_delay::Int64, integrator, rng)
    u = integrator.u
    @inbounds for (i, ξ) in delay_complete
        u[i] += num_next_delay*ξ
    end
end

function execute_delay_complete!(delay_complete::Function, num_next_delay::Int64, integrator, rng)
    for _ in 1:num_next_delay
        delay_complete(integrator, rng)
    end
end 


# """
#     function find_next_delay(de_chan::Vector{Vector{T}})

# Find the minimum
# """
# function find_next_delay(de_chan::Vector{Vector{T}}) where {T}
#     min_vec = Vector{T}(undef,length(de_chan))
#     for i in eachindex(de_chan)
#         min_vec[i] = isempty(de_chan[i]) ? Inf : minimum(de_chan[i])
#     end
#     argmin(min_vec)
# end

# """
# check number of next delay
# """
# function check_num_next_delay(one_de_chan,dt_delay)
#     count(i->(i==dt_delay),one_de_chan)
# end


# """
# find the minimal dt_delay in various delay channel
# """

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
function find_next_delay_vec(A::Vector{Vector{T}}, x::T) where {T}
    position_index = findall(A->x in A, A)
    num_in_vec = find_num_in_vec(A, position_index, x)
    return position_index, num_in_vec
end
