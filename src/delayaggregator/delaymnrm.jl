mutable struct DelayMNRMJumpAggregation{T,S,F1,F2,RNG,DG,PQ} <: AbstractDSSAJumpAggregator
    next_jump::Int
    prev_jump::Int
    next_jump_time::T
    end_time::T
    cur_rates::Vector{T}
    sum_rate::T
    ma_jumps::S
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
    rng::RNG
    dep_gr::DG
    pq::PQ
    next_delay::Union{Nothing,Vector{Int}}
    num_next_delay::Union{Nothing,Vector{Int}}
    time_to_next_jump::T
    dt_delay::T
end

function DelayMNRMJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; num_specs, dep_graph=nothing, kwargs...) where {T,S,F1,F2,RNG}
    

    # a dependency graph is needed and must be provided if there are constant rate jumps
    if dep_graph === nothing
        if (get_num_majumps(maj) == 0) || !isempty(rs)
            error("To use ConstantRateJumps with the Next Reaction Method (NRM) algorithm a dependency graph must be supplied.")
        else
            dg = make_dependency_graph(num_specs, maj)
        end
    else
        dg = dep_graph
        # make sure each jump depends on itself
        add_self_dependencies!(dg)
    end

    pq = MutableBinaryMinHeap{T}()

    nd = nothing
    nnd = nothing
    ttnj = zero(et)
    dt_delay = zero(et)
    DelayMNRMJumpAggregation{T,S,F1,F2,RNG,typeof(dg),typeof(pq)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, dg, pq, nd, nnd, ttnj, dt_delay)
end

############################# Required Functions ##############################
# creating the JumpAggregation structure (function wrapper-based constant jumps)
function aggregate(aggregator::DelayMNRM, u, p, t, end_time, constant_jumps,
    ma_jumps, save_positions, rng; kwargs...)

# handle constant jumps using function wrappers
rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)

build_jump_aggregation(DelayMNRMJumpAggregation, u, p, t, end_time, ma_jumps,
            rates, affects!, save_positions, rng; num_specs=length(u), kwargs...)
end



# set up a new simulation and calculate the first jump / jump time
function initialize!(p::DelayMNRMJumpAggregation, integrator, u, params, t)
    fill_rates_and_get_times!(p, u, params, t)
    
    find_next_delay_dt!(p, integrator)
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

# execute one jump, changing the system state
function execute_jumps!(p::DelayMNRMJumpAggregation, integrator, u, params, t)
    # execute jump
    update_state_delay!(p, integrator, u, t)
    # update current jump rates and times
    update_dependent_rates_delay!(p, integrator, u, params, t)
    nothing
end

# calculate the next jump / jump time
# just the top of the priority queue
function generate_jumps!(p::DelayMNRMJumpAggregation, integrator, u, params, t)
    next_jump_time_reaction, next_jump = top_with_handle(p.pq)
    dt_reaction = next_jump_time_reaction - t
    dt_delay_generation!(p, integrator)
    compare_delay!(p, integrator.de_chan, p.dt_delay, dt_reaction, t)
    if p.next_delay != nothing
        p.next_jump =0
    else
        p.next_jump = next_jump
    end
    nothing
end


# recalculate jump rates for jumps that depend on the just executed jump (p.next_jump)
function update_dependent_rates_delay!(p::DelayMNRMJumpAggregation, integrator, u, params, t)
    if p.next_delay == nothing  # if next reaction is not delay reaction 
        @inbounds dep_rxs = p.dep_gr[p.next_jump]
    else
        # find the dep_rxs w.r.t next_delay vectors
        vars_ = Vector{Vector{Int}}(undef,length(p.next_delay))
        var_to_jumps = var_to_jumps_map(length(u),p.ma_jumps)
        @inbounds for i in eachindex(p.next_delay)
            vars_[i] = first.(integrator.delayjumpsets.delay_complete[p.next_delay[i]])
        end
        vars = reduce(vcat, vars_)
        dep_rxs_ = Vector{Vector{Int}}(undef,length(vars))
        @inbounds for i in eachindex(dep_rxs_)
            dep_rxs_[i] = var_to_jumps[vars[i]]
        end
        dep_rxs = reduce(vcat, dep_rxs_)
    end
    @unpack cur_rates, rates, ma_jumps = p
    num_majumps = get_num_majumps(ma_jumps)
    @inbounds for rx in dep_rxs
        oldrate = cur_rates[rx]

        # update the jump rate
        @inbounds cur_rates[rx] = calculate_jump_rate(ma_jumps, num_majumps, rates, u, params, t, rx)

        # calculate new jump times for dependent jumps
        if rx != p.next_jump && oldrate > zero(oldrate)
            if cur_rates[rx] > zero(eltype(cur_rates))
                update!(p.pq, rx, t + oldrate / cur_rates[rx] * (p.pq[rx] - t))
            else
                update!(p.pq, rx, typemax(t))
            end
        else
            if cur_rates[rx] > zero(eltype(cur_rates))
                update!(p.pq, rx, t + randexp(p.rng) / cur_rates[rx])
            else
                update!(p.pq, rx, typemax(t))
            end
        end
    end
    nothing
end

# reevaulate all rates, recalculate all jump times, and reinit the priority queue
function fill_rates_and_get_times!(p::DelayMNRMJumpAggregation, u, params, t)

    # mass action jumps
    majumps   = p.ma_jumps
    cur_rates = p.cur_rates
    # num_de_chan = length(integrator.de_chan) # num_de_chan
    pqdata = Vector{typeof(t)}(undef,length(cur_rates))
    @inbounds for i in 1:get_num_majumps(majumps)
        cur_rates[i] = evalrxrate(u, i, majumps)
        pqdata[i] = t + randexp(p.rng) / cur_rates[i]
    end

    # constant rates
    rates = p.rates
    idx   = get_num_majumps(majumps) + 1
    @inbounds for rate in rates
        cur_rates[idx] = rate(u, params, t)
        pqdata[idx] = t + randexp(p.rng) / cur_rates[idx]
        idx += 1
    end

    # setup a new indexed priority queue to storing rx times
    p.pq = MutableBinaryMinHeap(pqdata)
    nothing
end
