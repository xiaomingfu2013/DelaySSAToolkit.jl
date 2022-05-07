mutable struct DelayDirectJumpAggregation{T,S,F1,F2,RNG,IType} <: AbstractDSSAJumpAggregator
    next_jump::Int
    prev_jump::Int
    next_jump_time::T
    end_time::T
    cur_rates::Vector{T} # cur_rates here is the cumsum of rates
    sum_rate::T
    ma_jumps::S
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
    rng::RNG
    next_delay::Union{Nothing,Vector{Int}}
    num_next_delay::Union{Nothing,Vector{Int}}
    time_to_next_jump::T
    shadow_integrator::IType
end
mutable struct ShadowIntegrator{uType,chanType,T}
    u::uType
    de_chan::chanType
    delayjumpsets::DelayJumpSet
    cur_rates::Vector{T}
    sum_rate::T
end
function DelayDirectJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; u0, kwargs...) where {T,S,F1,F2,RNG}
    ttnj = zero(et)
    shadow_integrator = ShadowIntegrator{typeof(u0),Vector{Vector{T}},T}(copy(u0), [Vector{T}()], DelayJumpSet(Dict(), Dict(), Dict()), copy(crs), sr)
    nd = nothing
    nnd = [1] # in the Direct Method the number of next delay equals always 1
    DelayDirectJumpAggregation{T,S,F1,F2,RNG,typeof(shadow_integrator)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, nd, nnd, ttnj, shadow_integrator)
end

function aggregate(aggregator::DelayDirect, u, p, t, end_time, constant_jumps, ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using function wrappers
    rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)
    build_jump_aggregation(DelayDirectJumpAggregation, u, p, t, end_time, ma_jumps, rates, affects!, save_positions, rng; u0=copy(u), kwargs...)
end

function initialize!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    p.shadow_integrator.delayjumpsets = integrator.delayjumpsets
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

@inline function execute_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    update_state_delay_direct!(p, integrator, integrator.u, t)
    nothing
end

function generate_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    generate_time_to_next_jump!(p, integrator, params, t)
    @fastmath p.next_jump_time = t + p.time_to_next_jump
    @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rand(p.rng) * p.sum_rate) # 
    nothing
end
"""
    Create delta based on the shawdow variable u_shadow
"""
@inline function generate_time_to_next_jump!(p::DelayDirectJumpAggregation, integrator, params, t)
    # the reason to use a shadow_integrator is because generating ttnj and changing the state happen simultaneously in DelayDirect method, so has to cache it before execute_jumps!
    p.shadow_integrator.de_chan = deepcopy(integrator.de_chan) #TODO
    p.shadow_integrator.u = copy(integrator.u) #TODO

    direct_algo!(p, p.shadow_integrator, params, t)
end

function direct_algo!(p, shadow_integrator, params, t; tgap = typemax(t))
    calculate_sum_rate!(p, shadow_integrator, shadow_integrator.u, params, t)
    r1 = rand(p.rng)
    if isempty(reduce(vcat, shadow_integrator.de_chan))
        ttnj = -log(r1) / shadow_integrator.sum_rate
        ttnj_last = ttnj
    else
        prev_T1 = zero(t)
        cur_T1, cur_T2 = find_next_delay_num(shadow_integrator.de_chan)
        i = 1
        aₜ = shadow_integrator.sum_rate * cur_T1
        F = one(t) - exp(-aₜ)
        aₜ_ = zero(aₜ)
        # shadow_integrator.sum_rate = p.sum_rate
        while F < r1 && prev_T1 <= tgap
            p.next_delay = [cur_T2]

            shift_delay_channel!(shadow_integrator.de_chan, cur_T1 - prev_T1)
            update_delay_channel!(shadow_integrator.de_chan)
            update_delay_complete!(p, shadow_integrator)

            # add support to handle T that is changing
            calculate_sum_rate!(p, shadow_integrator, shadow_integrator.u, params, t + cur_T1)

            prev_T1 = cur_T1
            cur_T1_, cur_T2 = find_next_delay_num(shadow_integrator.de_chan) 
            cur_T1 = cur_T1_ + prev_T1

            # aₜ_ = copy(aₜ) # backup aₜ
            aₜ_ = aₜ # backup aₜ
            aₜ += shadow_integrator.sum_rate * (cur_T1 - prev_T1)
            F = one(t) - exp(-aₜ)
            i += 1
        end
        ttnj_last = tgap<typemax(t) ? tgap - prev_T1 : (-log(one(t) - r1) - aₜ_) / shadow_integrator.sum_rate
        ttnj = prev_T1 + ttnj_last
    end
    
    # T1_last, T2_last = create_Tstruct(integrator.de_chan)

    # ttnj_last will not change the state anymore
    shift_delay_channel!(shadow_integrator.de_chan, ttnj_last)
    update_delay_channel!(shadow_integrator.de_chan)

    # in case the last ttnj also change the state
    # update_state_final_jump!(p, integrator, ttnj_last, T1_last, T2_last)
    if tgap == typemax(t)
        fill_cum_rates_and_sum!(p, shadow_integrator.u, params, t + ttnj)
        p.time_to_next_jump = ttnj
    end
end

# @inbounds function update_state_final_jump!(p, integrator, tgap, T1, T2)
#     idx = count(x -> x <= tgap, T1)
#     for i in 1:idx
#         p.next_delay = [T2[i]]
#         update_delay_complete!(p, integrator)
#     end
#     nothing
# end

function update_delay_chan_state_at_tstop_test!(p, integrator, params, t, tgap)
    p.shadow_integrator.u = integrator.u
    p.shadow_integrator.de_chan = integrator.de_chan
    direct_algo!(p, p.shadow_integrator, params, t; tgap = tgap)
    nothing
end

"""
    function fill_cum_rates_and_sum!(p::DelayDirectJumpAggregation, u, params, t)

This function modifies `p.cur_rates` and `p.sum_rate`
"""
function fill_cum_rates_and_sum!(p::DelayDirectJumpAggregation, u, params, t)
    prev_rate = zero(t)
    new_rate = zero(t)
    cur_rates = p.cur_rates
    # mass action rates
    majumps = p.ma_jumps
    idx = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate = cur_rates[i]
    end
    # constant jump rates
    idx += 1
    rates = p.rates
    @inbounds for i in 1:length(p.rates)
        new_rate = rates[i](u, params, t)
        cur_rates[idx] = new_rate + prev_rate
        prev_rate = cur_rates[idx]
        idx += 1
    end

    @inbounds sum_rate = cur_rates[end]
    p.sum_rate = sum_rate
    nothing
end
"""
    Only changes s.sum_rate
"""
function calculate_sum_rate!(p, s::ShadowIntegrator, u, params, t)
    prev_rate = zero(t)
    new_rate = zero(t)
    cur_rates = s.cur_rates
    majumps = p.ma_jumps
    idx = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate = cur_rates[i]
    end
    # constant jump rates
    idx += 1
    rates = p.rates
    @inbounds for i in 1:length(rates)
        new_rate = rates[i](u, params, t)
        cur_rates[idx] = new_rate + prev_rate
        prev_rate = cur_rates[idx]
        idx += 1
    end
    s.sum_rate = cur_rates[end]
end

@inline function update_state_delay_direct!(p::DelayDirectJumpAggregation, integrator, u, t)
    @unpack ma_jumps, next_jump, time_to_next_jump = p
    @unpack delay_trigger_set, delay_interrupt_set = integrator.delayjumpsets

    integrator.u = copy(p.shadow_integrator.u)
    integrator.de_chan = deepcopy(p.shadow_integrator.de_chan) #TODO 

    num_ma_rates = get_num_majumps(ma_jumps)
    if next_jump <= num_ma_rates # if the next jump is a mass action jump
        if u isa SVector
            integrator.u = executerx(u, next_jump, ma_jumps)
        else
            @inbounds executerx!(integrator.u, next_jump, ma_jumps)
        end
    else
        idx = next_jump - num_ma_rates
        @inbounds p.affects![idx](integrator)
    end

    if next_jump in delay_interrupt_set
        update_delay_interrupt!(p, integrator)
    elseif next_jump in delay_trigger_set
        update_delay_trigger!(p, integrator)
    end
    # save jump that was just executed
    p.prev_jump = next_jump
    nothing
end

