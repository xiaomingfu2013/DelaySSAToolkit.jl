const MINJUMPRATE = 2.0^exponent(1e-12)

mutable struct DelayDirectCRJumpAggregation{T,S,F1,F2,RNG,DEPGR,U<:DiffEqJump.PriorityTable,W<:Function} <: AbstractDSSAJumpAggregator
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
    dep_gr::DEPGR
    minrate::T
    maxrate::T   # initial maxrate only, table can increase beyond it!
    rt::U
    ratetogroup::W
    next_delay::Vector{Int}
    num_next_delay::Vector{Int}
    time_to_next_jump::T
    dt_delay::T
    vartojumps_map::Union{Nothing,Vector{Vector{Int}}}
    dep_gr_delay::Union{Nothing,Dict{Int,Vector{Int}}}
end

function DelayDirectCRJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T,
    maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool},
    rng::RNG; num_specs, dep_graph=nothing, dep_graph_delay=nothing,
    minrate=convert(T, MINJUMPRATE), maxrate=convert(T, Inf), vartojumps_map = nothing,
    kwargs...) where {T,S,F1,F2,RNG}

    # a dependency graph is needed and must be provided if there are constant rate jumps
    if dep_graph === nothing
        if (get_num_majumps(maj) == 0) || !isempty(rs)
            error("To use ConstantRateJumps with the DelayDirectCR algorithm a dependency graph must be supplied.")
        else
            dg = make_dependency_graph(num_specs, maj)
        end
    else
        dg = dep_graph

        # make sure each jump depends on itself
        add_self_dependencies!(dg)
    end

    # mapping from jump rate to group id
    minexponent = exponent(minrate)

    # use the largest power of two that is <= the passed in minrate
    minrate = 2.0^minexponent
    ratetogroup = rate -> DiffEqJump.priortogid(rate, minexponent)

    # construct an empty initial priority table -- we'll reset this in init
    rt = DiffEqJump.PriorityTable(ratetogroup, zeros(T, 1), minrate, 2 * minrate)
    nd = Int64[]
    nnd = Int64[]
    ttnj = zero(et)
    dt_delay = zero(et)
    if vartojumps_map === nothing
        if (get_num_majumps(maj) == 0) || !isempty(rs)
            if dep_graph_delay === nothing
                @warn "To use ConstantRateJumps with the DelayDirectCR algorithm: make sure a delay dependency graph is correctly supplied!"
                vartojumps_map = repeat([1:length(crs)], num_specs)
            end
        else
            vartojumps_map = var_to_jumps_map(num_specs, maj)
        end
    end
    dep_gr_delay = dep_graph_delay
    DelayDirectCRJumpAggregation{T,S,F1,F2,RNG,typeof(dg),typeof(rt),typeof(ratetogroup)}(
        nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng,
        dg, minrate, maxrate, rt, ratetogroup, nd, nnd, ttnj, dt_delay, vartojumps_map, dep_gr_delay)
end


############################# Required Functions ##############################

# creating the JumpAggregation structure (function wrapper-based constant jumps)
function aggregate(aggregator::DelayDirectCR, u, p, t, end_time, constant_jumps,
    ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using function wrappers
    rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)

    build_jump_aggregation(DelayDirectCRJumpAggregation, u, p, t, end_time, ma_jumps,
        rates, affects!, save_positions, rng; num_specs=length(u), kwargs...)
end


# set up a new simulation and calculate the first jump / jump time
function initialize!(p::DelayDirectCRJumpAggregation, integrator, u, params, t)

    # initialize rates
    fill_rates_and_sum!(p, u, params, t)
    if p.dep_gr_delay === nothing
        p.dep_gr_delay = dep_gr_delay(integrator.delayjumpsets, p.vartojumps_map, length(p.cur_rates))
    end
    # setup PriorityTable
    DiffEqJump.reset!(p.rt)
    for (pid, priority) in enumerate(p.cur_rates)
        DiffEqJump.insert!(p.rt, pid, priority)
    end
    find_next_delay_dt!(p, integrator)
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

# execute one jump, changing the system state
function execute_jumps!(p::DelayDirectCRJumpAggregation, integrator, u, params, t)
    # execute jump
    update_state_delay!(p, integrator, u, t)

    # update current jump rates
    update_dependent_rates_delay!(p, integrator, integrator.u, params, t)
    nothing
end

# calculate the next jump / jump time
function generate_jumps!(p::DelayDirectCRJumpAggregation, integrator, u, params, t)
    dt_reaction = randexp(p.rng) / p.sum_rate
    dt_delay_generation!(p, integrator)
    compare_delay!(p, integrator.de_chan, p.dt_delay, dt_reaction, t)
    if p.next_jump_time < p.end_time
        if !isempty(p.next_delay)
            p.next_jump = 0
        else
            p.next_jump = DiffEqJump.sample(p.rt, p.cur_rates, p.rng)
        end
    end
    nothing
end




######################## SSA specific helper routines #########################

# recalculate jump rates for jumps that depend on the just executed jump
# requires dependency graph
function update_dependent_rates_delay!(p::DelayDirectCRJumpAggregation, integrator, u, params, t)

    if isempty(p.next_delay)  # if next reaction is not delay reaction 
        @inbounds dep_rxs = p.dep_gr[p.next_jump]
    else
        # find the dep_rxs w.r.t next_delay vectors
        dep_rxs_ = [p.dep_gr_delay[p.next_delay[i]] for i in eachindex(p.next_delay)]
        dep_rxs = reduce(vcat, dep_rxs_)
    end

    @unpack cur_rates, rates, ma_jumps, rt = p
    num_majumps = get_num_majumps(ma_jumps)

    @inbounds for rx in dep_rxs
        oldrate = cur_rates[rx]

        # update rate
        cur_rates[rx] = calculate_jump_rate(ma_jumps, num_majumps, rates, u, params, t, rx)

        # update table
        DiffEqJump.update!(rt, rx, oldrate, cur_rates[rx])
    end

    p.sum_rate = DiffEqJump.groupsum(rt)
    nothing
end