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
    next_delay::Vector{Int}
    num_next_delay::Vector{Int}
    time_to_next_jump::T
    dt_delay::T
    vartojumps_map::Union{Nothing,Vector{Vector{Int}}}
    dep_gr_delay::Union{Nothing,Dict{Int,Vector{Int}}}
end

function DelayMNRMJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; num_specs, dep_graph = nothing, dep_graph_delay = nothing, vartojumps_map = nothing, kwargs...) where {T,S,F1,F2,RNG}


    # a dependency graph is needed and must be provided if there are constant rate jumps
    if dep_graph === nothing
        if (get_num_majumps(maj) == 0) || !isempty(rs)
            error("To use ConstantRateJumps with the Delay Modified Next Reaction Method (Delay MNRM) algorithm a dependency graph must be supplied.")
        else
            dg = make_dependency_graph(num_specs, maj)
        end
    else
        dg = dep_graph
        # make sure each jump depends on itself
        add_self_dependencies!(dg)
    end


    pq = MutableBinaryMinHeap{T}()

    nd = Int64[]
    nnd = Int64[]
    ttnj = zero(et)
    dt_delay = zero(et)
    if vartojumps_map === nothing
        if (get_num_majumps(maj) == 0) || !isempty(rs)
            if dep_graph_delay === nothing  
                @warn "To use ConstantRateJumps with the DelayMNRM algorithm: make sure a delay dependency graph is correctly supplied!"
                vartojumps_map = repeat([1:length(crs)], num_specs)
            end
        else
            vartojumps_map = var_to_jumps_map(num_specs, maj)
        end
    end
    dep_gr_delay = dep_graph_delay
    DelayMNRMJumpAggregation{T,S,F1,F2,RNG,typeof(dg),typeof(pq)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, dg, pq, nd, nnd, ttnj, dt_delay, vartojumps_map, dep_gr_delay)
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
    if p.dep_gr_delay === nothing
        p.dep_gr_delay = dep_gr_delay(integrator.delayjumpsets, p.vartojumps_map, length(p.cur_rates))
    end
    find_next_delay_dt!(p, integrator)
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

# execute one jump, changing the system state
function execute_jumps!(p::DelayMNRMJumpAggregation, integrator, u, params, t)
    # execute jump
    update_state_delay!(p, integrator, u, t)
    # update current jump rates and times
    update_dependent_rates_delay!(p, integrator, integrator.u, params, t)
    nothing
end

# calculate the next jump / jump time
# just the top of the priority queue
function generate_jumps!(p::DelayMNRMJumpAggregation, integrator, u, params, t)
    next_jump_time_reaction, next_jump = top_with_handle(p.pq)
    dt_reaction = next_jump_time_reaction - t
    dt_delay_generation!(p, integrator)
    compare_delay!(p, integrator.de_chan, p.dt_delay, dt_reaction, t)
    if !isempty(p.next_delay)
        p.next_jump = 0
    else
        p.next_jump = next_jump
    end
    nothing
end


# recalculate jump rates for jumps that depend on the just executed jump (p.next_jump)
function update_dependent_rates_delay!(p::DelayMNRMJumpAggregation, integrator, u, params, t)
    if isempty(p.next_delay)  # if next reaction is not delay reaction 
        @inbounds dep_rxs = p.dep_gr[p.next_jump]
    else
        # find the dep_rxs w.r.t next_delay vectors
        dep_rxs_ = [p.dep_gr_delay[p.next_delay[i]] for i in eachindex(p.next_delay)]
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
                DataStructures.update!(p.pq, rx, t + oldrate / cur_rates[rx] * (p.pq[rx] - t))
            else
                DataStructures.update!(p.pq, rx, typemax(t))
            end
        else
            if cur_rates[rx] > zero(eltype(cur_rates))
                DataStructures.update!(p.pq, rx, t + randexp(p.rng) / cur_rates[rx])
            else
                DataStructures.update!(p.pq, rx, typemax(t))
            end
        end
    end
    nothing
end

# reevaulate all rates, recalculate all jump times, and reinit the priority queue
function fill_rates_and_get_times!(p::DelayMNRMJumpAggregation, u, params, t)

    # mass action jumps
    majumps = p.ma_jumps
    cur_rates = p.cur_rates
    # num_de_chan = length(integrator.de_chan) # num_de_chan
    pqdata = Vector{typeof(t)}(undef, length(cur_rates))
    @inbounds for i in 1:get_num_majumps(majumps)
        cur_rates[i] = evalrxrate(u, i, majumps)
        pqdata[i] = t + randexp(p.rng) / cur_rates[i]
    end

    # constant rates
    rates = p.rates
    idx = get_num_majumps(majumps) + 1
    @inbounds for rate in rates
        cur_rates[idx] = rate(u, params, t)
        pqdata[idx] = t + randexp(p.rng) / cur_rates[idx]
        idx += 1
    end

    # setup a new indexed priority queue to storing rx times
    p.pq = MutableBinaryMinHeap(pqdata)
    nothing
end
