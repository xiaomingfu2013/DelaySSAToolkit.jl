"""
Queue method. This method handles variable intensity rates.
See https://github.com/SciML/JumpProcesses.jl/pull/276
"""
mutable struct DelayCoevolveJumpAggregation{T, S, F1, F2, RNG, GR, PQ} <:
               AbstractDSSAJumpAggregator
    next_jump::Int                    # the next jump to execute
    prev_jump::Int                    # the previous jump that was executed
    next_jump_time::T                 # the time of the next jump
    end_time::T                       # the time to stop a simulation
    cur_rates::Vector{T}              # the last computed upper bound for each rate
    sum_rate::Nothing                 # not used
    ma_jumps::S                       # MassActionJumps
    rates::F1                         # vector of rate functions
    affects!::F2                      # vector of affect functions for VariableRateJumps
    save_positions::Tuple{Bool, Bool} # tuple for whether to save the jumps before and/or after event
    rng::RNG                          # random number generator
    dep_gr::GR                        # map from jumps to jumps depending on it
    pq::PQ                            # priority queue of next time
    lrates::F1                        # vector of rate lower bound functions
    urates::F1                        # vector of rate upper bound functions
    rateintervals::F1                 # vector of interval length functions
    haslratevec::Vector{Bool}         # vector of whether an lrate was provided for this vrj
    next_delay::Vector{Int}
    num_next_delay::Vector{Int}
    time_to_next_jump::T
    dt_delay::T
    vartojumps_map::Union{Nothing, Vector{Vector{Int}}}
    dep_gr_delay::Union{Nothing, Dict{Int, Vector{Int}}}
end

function DelayCoevolveJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::Nothing,
                                      maj::S, rs::F1, affs!::F2, sps::Tuple{Bool, Bool},
                                      rng::RNG; u::U, dep_graph = nothing, lrates, urates,
                                      rateintervals, haslratevec, dep_graph_delay = nothing,
                                      vartojumps_map = nothing) where {T, S, F1, F2, RNG, U}
    if dep_graph === nothing
        if (get_num_majumps(maj) == 0) || !isempty(urates)
            error("To use DelayCoevolve a dependency graph between jumps must be supplied.")
        else
            dg = make_dependency_graph(length(u), maj)
        end
    else
        # using a Set to ensure that edges are not duplicate
        dgsets = [Set{Int}(append!(Int[], jumps, [var]))
                  for (var, jumps) in enumerate(dep_graph)]
        dg = [sort!(collect(i)) for i in dgsets]
    end

    num_jumps = get_num_majumps(maj) + length(urates)

    if length(dg) != num_jumps
        error("Number of nodes in the dependency graph must be the same as the number of jumps.")
    end

    pq = MutableBinaryMinHeap{T}()

    num_specs = length(u)

    nd = Int64[]
    nnd = Int64[]
    ttnj = zero(et)
    dt_delay = zero(et)
    if vartojumps_map === nothing
        # println(dep_graph_delay)
        # println(urates)
        if (get_num_majumps(maj) == 0) || !isempty(urates)
            if dep_graph_delay === nothing
                @warn "To use ConstantRateJumps with this algorithm: make sure a delay dependency graph is correctly supplied!"
                vartojumps_map = repeat([1:length(crs)], num_specs)
            end
        else
            vartojumps_map = var_to_jumps_map(num_specs, maj)
        end
    end
    dep_gr_delay = dep_graph_delay
    DelayCoevolveJumpAggregation{T, S, F1, F2, RNG, typeof(dg),
                                 typeof(pq)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps,
                                             rng,
                                             dg, pq, lrates, urates, rateintervals,
                                             haslratevec, nd, nnd, ttnj, dt_delay,
                                             vartojumps_map, dep_gr_delay)
end

# creating the JumpAggregation structure (tuple-based variable jumps)
function aggregate(aggregator::DelayCoevolve, u, p, t, end_time, constant_jumps,
                   ma_jumps, save_positions, rng; dep_graph = nothing,
                   variable_jumps = nothing, dep_graph_delay = nothing,
                   vartojumps_map = nothing, kwargs...)
    AffectWrapper = FunctionWrappers.FunctionWrapper{Nothing, Tuple{Any}}
    RateWrapper = FunctionWrappers.FunctionWrapper{typeof(t),
                                                   Tuple{typeof(u), typeof(p), typeof(t)}}

    ncrjs = (constant_jumps === nothing) ? 0 : length(constant_jumps)
    nvrjs = (variable_jumps === nothing) ? 0 : length(variable_jumps)
    nrjs = ncrjs + nvrjs
    affects! = Vector{AffectWrapper}(undef, nrjs)
    rates = Vector{RateWrapper}(undef, nvrjs)
    lrates = similar(rates)
    rateintervals = similar(rates)
    urates = Vector{RateWrapper}(undef, nrjs)
    haslratevec = zeros(Bool, nvrjs)

    idx = 1
    if constant_jumps !== nothing
        for crj in constant_jumps
            affects![idx] = AffectWrapper(integ -> (crj.affect!(integ); nothing))
            urates[idx] = RateWrapper(crj.rate)
            idx += 1
        end
    end

    if variable_jumps !== nothing
        for (i, vrj) in enumerate(variable_jumps)
            affects![idx] = AffectWrapper(integ -> (vrj.affect!(integ); nothing))
            urates[idx] = RateWrapper(vrj.urate)
            idx += 1
            rates[i] = RateWrapper(vrj.rate)
            rateintervals[i] = RateWrapper(vrj.rateinterval)
            haslratevec[i] = haslrate(vrj)
            lrates[i] = haslratevec[i] ? RateWrapper(vrj.lrate) : RateWrapper(nullrate)
        end
    end

    num_jumps = get_num_majumps(ma_jumps) + nrjs
    cur_rates = Vector{typeof(t)}(undef, num_jumps)
    sum_rate = nothing
    next_jump = 0
    next_jump_time = typemax(t)
    DelayCoevolveJumpAggregation(next_jump, next_jump_time, end_time, cur_rates, sum_rate,
                                 ma_jumps, rates, affects!, save_positions, rng;
                                 u, dep_graph, lrates, urates, rateintervals, haslratevec,
                                 dep_graph_delay, vartojumps_map)
end

# set up a new simulation and calculate the first jump / jump time
function initialize!(p::DelayCoevolveJumpAggregation, integrator, u, params, t)
    p.end_time = integrator.sol.prob.tspan[2]
    fill_rates_and_get_times!(p, u, params, t)
    if p.dep_gr_delay === nothing
        p.dep_gr_delay = dep_gr_delay(integrator.delayjumpsets, p.vartojumps_map,
                                      length(p.cur_rates))
    end
    find_next_delay_dt!(p, integrator)
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

# execute one jump, changing the system state
function execute_jumps!(p::DelayCoevolveJumpAggregation, integrator, u, params, t)
    # execute jump
    update_state_delay!(p, integrator, u, t)
    # update current jump rates and times
    # update_dependent_rates!(p, u, params, t)
    update_dependent_rates_delay!(p, integrator, integrator.u, params, t)
    nothing
end

# calculate the next jump / jump time
function generate_jumps!(p::DelayCoevolveJumpAggregation, integrator, u, params, t)
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

######################## SSA specific helper routines ########################
function update_dependent_rates_delay!(p::DelayCoevolveJumpAggregation, integrator, u,
                                       params, t)
    if isempty(p.next_delay)  # if next reaction is not delay reaction 
        @inbounds deps = p.dep_gr[p.next_jump]
    else
        # find the dep_rxs w.r.t next_delay vectors
        dep_rxs_ = [p.dep_gr_delay[p.next_delay[i]] for i in eachindex(p.next_delay)]
        deps = reduce(vcat, dep_rxs_)
    end
    @unpack cur_rates, end_time, pq = p
    for (ix, i) in enumerate(deps)
        ti, last_urate_i = next_time(p, u, params, t, i, end_time)
        update!(pq, i, ti)
        @inbounds cur_rates[i] = last_urate_i
    end
    nothing
end

@inline function get_ma_urate(p::DelayCoevolveJumpAggregation, i, u, params, t)
    return evalrxrate(u, i, p.ma_jumps)
end

@inline function get_urate(p::DelayCoevolveJumpAggregation, uidx, u, params, t)
    @inbounds return p.urates[uidx](u, params, t)
end

@inline function get_rateinterval(p::DelayCoevolveJumpAggregation, lidx, u, params, t)
    @inbounds return p.rateintervals[lidx](u, params, t)
end

@inline function get_lrate(p::DelayCoevolveJumpAggregation, lidx, u, params, t)
    @inbounds return p.lrates[lidx](u, params, t)
end

@inline function get_rate(p::DelayCoevolveJumpAggregation, lidx, u, params, t)
    @inbounds return p.rates[lidx](u, params, t)
end

function next_time(p::DelayCoevolveJumpAggregation{T}, u, params, t, i, tstop::T) where {T}
    @unpack rng, haslratevec = p
    num_majumps = get_num_majumps(p.ma_jumps)
    num_cjumps = length(p.urates) - length(p.rates)
    uidx = i - num_majumps
    lidx = uidx - num_cjumps
    urate = uidx > 0 ? get_urate(p, uidx, u, params, t) : get_ma_urate(p, i, u, params, t)
    last_urate = p.cur_rates[i]
    if i != p.next_jump && last_urate > zero(t)
        s = urate == zero(t) ? typemax(t) : last_urate / urate * (p.pq[i] - t)
    else
        s = urate == zero(t) ? typemax(t) : randexp(rng) / urate
    end
    _t = t + s
    if lidx > 0
        while t < tstop
            rateinterval = get_rateinterval(p, lidx, u, params, t)
            if s > rateinterval
                t = t + rateinterval
                urate = get_urate(p, uidx, u, params, t)
                s = urate == zero(t) ? typemax(t) : randexp(rng) / urate
                _t = t + s
                continue
            end
            (_t >= tstop) && break

            lrate = haslratevec[lidx] ? get_lrate(p, lidx, u, params, t) : zero(t)
            if lrate < urate
                # when the lower and upper bound are the same, then v < 1 = lrate / urate = urate / urate
                v = rand(rng) * urate
                # first inequality is less expensive and short-circuits the evaluation
                if (v > lrate) && (v > get_rate(p, lidx, u, params, _t))
                    t = _t
                    urate = get_urate(p, uidx, u, params, t)
                    s = urate == zero(t) ? typemax(t) : randexp(rng) / urate
                    _t = t + s
                    continue
                end
            elseif lrate > urate
                error("The lower bound should be lower than the upper bound rate for t = $(t) and i = $(i), but lower bound = $(lrate) > upper bound = $(urate)")
            end
            break
        end
    end
    return _t, urate
end

# reevaulate all rates, recalculate all jump times, and reinit the priority queue
function fill_rates_and_get_times!(p::DelayCoevolveJumpAggregation, u, params, t)
    @unpack end_time = p
    num_jumps = get_num_majumps(p.ma_jumps) + length(p.urates)
    p.cur_rates = zeros(typeof(t), num_jumps)
    jump_times = Vector{typeof(t)}(undef, num_jumps)
    @inbounds for i in 1:num_jumps
        jump_times[i], p.cur_rates[i] = next_time(p, u, params, t, i, end_time)
    end
    p.pq = MutableBinaryMinHeap(jump_times)
    nothing
end
