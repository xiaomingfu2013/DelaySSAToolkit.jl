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
    build_jump_aggregation(jump_agg_type, u, p, t, end_time, ma_jumps, rates,
                           affects!, save_positions, rng; kwargs...)

Helper routine for setting up standard fields of DSSA jump aggregations.
"""
function build_jump_aggregation(jump_agg_type, u, p, t, end_time, ma_jumps, rates, affects!, save_positions, rng; kwargs...)

    # mass action jumps
    majumps = ma_jumps
    if majumps === nothing
        majumps = MassActionJump(Vector{typeof(t)}(),Vector{Vector{Pair{Int,eltype(u)}}}(),Vector{Vector{Pair{Int,eltype(u)}}}())
    end

    # current jump rates, allows mass action rates and constant jumps
    cur_rates = Vector{typeof(t)}(undef, get_num_majumps(majumps) + length(rates))

    sum_rate = zero(typeof(t))
    next_jump = 0
    next_jump_time = typemax(typeof(t))
    jump_agg_type(next_jump, next_jump_time, end_time, cur_rates, sum_rate,
                majumps, rates, affects!, save_positions, rng; kwargs...)
end


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


@inline @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rxidx::S,majump::MassActionJump{U,V,W,X})::R where {T,S,R,U <: AbstractVector{R},V,W,X}
    val = one(T)
    @inbounds for specstoch in majump.reactant_stoch[rxidx]
        specpop = speciesvec[specstoch[1]]
        val    *= specpop
        @inbounds for _ = 2:specstoch[2]
            specpop -= one(specpop)
            val     *= specpop
        end
    end
    @inbounds return val * majump.scaled_rates[rxidx]
end

@inline @fastmath function executerx!(speciesvec::AbstractVector{T}, rxidx::S,majump::M) where {T,S,M <: DiffEqJump.AbstractMassActionJump}
    @inbounds net_stoch = majump.net_stoch[rxidx]
    @inbounds for specstoch in net_stoch
        speciesvec[specstoch[1]] += specstoch[2]
    end
    nothing
end

@inline @fastmath function executerx(speciesvec::SVector{T}, rxidx::S,majump::M) where {T,S,M <: DiffEqJump.AbstractMassActionJump}
    @inbounds net_stoch = majump.net_stoch[rxidx]
    @inbounds for specstoch in net_stoch
        speciesvec = setindex(speciesvec,speciesvec[specstoch[1]]+specstoch[2],specstoch[1])
    end
    speciesvec
end

@inline function dt_delay_generation!(p, integrator)
    next_jump = p.next_jump
    ttnj = p.time_to_next_jump
    @unpack delay_trigger_set, delay_interrupt_set = integrator.delayjumpsets
    if !isempty(p.next_delay) || any(set->next_jump in set, [delay_interrupt_set, delay_trigger_set])
        find_next_delay_dt!(p, integrator)
    else
        p.dt_delay -= ttnj
    end
end

@inline function update_state_delay!(p::AbstractDSSAJumpAggregator, integrator, u, t)
    @unpack ma_jumps, next_jump, next_delay, time_to_next_jump = p
    @unpack delay_interrupt, delay_trigger, delay_trigger_set, delay_interrupt_set, delay_complete = integrator.delayjumpsets

    ttnj = time_to_next_jump
    if isempty(next_delay)
        # update_state ! 
        num_ma_rates = get_num_majumps(ma_jumps)
        if next_jump <= num_ma_rates # is next jump a mass action jump
            if u isa SVector
                integrator.u = executerx(integrator.u, next_jump, ma_jumps)
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
        next_delay = Int64[]
        num_next_delay = Int64[]
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
    p.dt_delay = find_minimun_dt_delay(de_chan)
    nothing
end

function find_minimun_dt_delay(de_chan::Vector{Vector{T}}) where {T}
    val_vec = Vector{T}(undef,length(de_chan))
    @inbounds for i in eachindex(de_chan)
        val_vec[i] = isempty(de_chan[i]) ? typemax(T) : minimum(de_chan[i])
    end
    minimum(val_vec)
end

function find_next_delay_num!(p, de_chan::Vector{Vector{T}}) where {T}
    @label restart
    val = find_minimun_dt_delay(de_chan)
    if val < eps(T)
        shift_delay_channel!(de_chan, eps(T))
        update_delay_channel!(de_chan)
        @goto restart
    end
    p.next_delay, p.num_next_delay = find_next_delay_vec(de_chan, val)
    p.next_delay_time = val
    nothing
end

@inline function shift_delay_channel!(de_chan::Vector{Vector{T1}},ttnj::T2) where {T1<:Real,T2<:Real}
    @inbounds for idx in eachindex(de_chan)
        for j in eachindex(de_chan[idx])
            de_chan[idx][j] -= ttnj
        end
    end
end

@inline function update_delay_channel!(de_chan::Vector{Vector{T}}) where {T<:Real}
    @inbounds for idx in eachindex(de_chan)
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
            integrator.u = setindex(integrator.u, integrator.u[i] + num_next_delay*ξ, i)
            # integrator.u = setindex!!(integrator.u, integrator.u[i] + num_next_delay*ξ, i)
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
    function vars_to_jumps_delay(p, vars)
input::Vector{Int} species indices
output::Vector{Int} reactions need to be updated
"""

function vars_to_jumps_delay(vartojumps_map::Vector{Vector{Int}}, vars::Vector{Int})
    jumps = []
    for var in vars
        push!(jumps, vartojumps_map[var])
    end
    return reduce(vcat, jumps)
end

"""
    function dep_gr_delay(p, integrator)
    Generate delay dependency graph
    output::Dict next_delay idx => reactions need to be updated
"""
function dep_gr_delay(delayjumpsets::DelayJumpSet, vartojumps_map, num_reactions)
    dict_ = Dict{Int,Vector{Int}}()
    dict_complete = delayjumpsets.delay_complete
    @inbounds for key in keys(dict_complete)
        delay_complete_action = dict_complete[key]
        if typeof(delay_complete_action)<:Vector{Pair{Int,Int}}
            vars = first.(delay_complete_action)
            jumps = unique(vars_to_jumps_delay(vartojumps_map, vars))
        else
            jumps = vec(1:num_reactions)
        end
        push!(dict_, key=>jumps)
    end
    dict_
end

"""
    find_next_delay_vec(A::Vector{Vector{T}}, x::T)
Returns two vectors. The first is the position indices; `position_indices` of vector `A`, and the second is the vector `num_in_vec` composed of the occurrence of value `x`.
Used in finding the minimal dt_delay in various delay channel.
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
    position_indices = Vector{Int64}()
    num_in_vec = Vector{Int64}()
    for idx in eachindex(de_chan)
        if ttnj in de_chan[idx]
            append!(position_indices, idx)
            append!(num_in_vec, count(==(ttnj), de_chan[idx]))
        end
    end
    return position_indices, num_in_vec
end


## decrepit
# """
#     function find_num_in_vec(A::Vector{Vector{T}}, position_indices::Vector{Int64}, x::T)

# Find the occurrence of value `x` in each vector, according to the corresponding position indices specified by the element in the `position_indices` vector in the given vector `A`.

# # Examples
# ```julia-repl
# julia> A =  [[0.09,0.09,0.1],[0.3,0.09,0.1],[0.09]]
# 3-element Vector{Vector{Float64}}:
#  [0.09, 0.09, 0.1]
#  [0.3, 0.09, 0.1]
#  [0.09]
# julia> position_indices =  [1,2]
# 2-element Vector{Int64}:
#  1
#  2
# julia> find_num_in_vec(A, position_indices, 0.09)
# 2-element Vector{Int64}:
# 2
# 1
# ```
# """
# @inline function find_num_in_vec(A::Vector{Vector{T}}, position_indices::Vector{Int64}, x::T) where {T}
#     number_in_vec = Vector{Int64}(undef, length(position_indices))
#     @inbounds for i in eachindex(position_indices)
#         number_in_vec[i] = count(==(x),A[position_indices[i]])
#     end
#     return number_in_vec
# end
# function find_next_delay_vec(de_chan::Vector{Vector{T}}, ttnj::T) where {T}
#     position_indices = findall(de_chan->ttnj in de_chan, de_chan)
#     num_in_vec = find_num_in_vec(de_chan, position_indices, ttnj)
#     return position_indices, num_in_vec
# end