"""
#TODO

"""

mutable struct DelayDirectJumpAggregation{T,S,F1,F2,RNG,uType} <: AbstractDSSAJumpAggregator
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
    time_to_next_jump::T
    u_shadow::uType
end
function DelayDirectJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; u0, kwargs...) where {T,S,F1,F2,RNG}
    ttnj = zero(et)
    u_shadow = u0
    DelayDirectJumpAggregation{T,S,F1,F2,RNG,typeof(u0)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, ttnj, u_shadow)
end

function aggregate(aggregator::DelayDirect, u, p, t, end_time, constant_jumps,
    ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using function wrappers
    rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)
    build_jump_aggregation(DelayDirectJumpAggregation, u, p, t, end_time, ma_jumps,
                           rates, affects!, save_positions, rng; u0 = copy(u), kwargs...)
end

function initialize!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

@inline function execute_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    integrator.u = copy(p.u_shadow) # Update the delay complete reactions
    update_state_delay_Direct!(p, integrator, integrator.u, t)
    nothing
end

function generate_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    p.u_shadow = copy(u) #TODO
    generate_delta!(p, integrator, params, t)

    @fastmath p.next_jump_time = t + p.time_to_next_jump
    @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rand(p.rng) * p.sum_rate) # 
    nothing
end
"""
Create delta based on the shawdow variable u_shadow
"""
@fastmath function generate_delta!(p::DelayDirectJumpAggregation, integrator, params, t)
    @unpack delay_complete = integrator.delayjumpsets
    
    fill_cum_rates_and_sum!(p, p.u_shadow, params, t)
    r1 = rand(p.rng)
    T1, T2 = create_Tstruct(integrator.de_chan)
    if isempty(T1)
        ttnj = -log(r1)/p.sum_rate
    else
        prepend!(T1,0)
        append!(T1,Inf)
        i = 1
        aₜ = p.sum_rate*T1[2]
        F = 1 - exp(-aₜ)
        while F < r1
            update_delay_complete!(p.u_shadow, [T2[i]], [1], delay_complete)
            if i < length(T1)-2
                u_ = copy(p.u_shadow)
                update_delay_complete!(u_, [T2[i+1]], [1], delay_complete)
                sum_rate_ = calculate_sum_rate(p, u_,  params, t+T1[i+2])
                aₜ += sum_rate_*(T1[i+2]-T1[i+1])
                F = 1 - exp(-aₜ)
                i += 1
            else
                F = 1
            end
        end
        sum_rate = calculate_sum_rate(p, p.u_shadow,  params, t+T1[i])
        ttnj = T1[i]-(log(1-r1)+aₜ-sum_rate*(T1[i+1]-T1[i]))/sum_rate
    end
    fill_cum_rates_and_sum!(p, p.u_shadow, params, t+ttnj)
    p.time_to_next_jump = ttnj
end

function fill_cum_rates_and_sum!(p::DelayDirectJumpAggregation, u, params, t)
    prev_rate = zero(t)
    new_rate  = zero(t)
    cur_rates = p.cur_rates
    # mass action rates
    majumps   = p.ma_jumps
    idx       = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate     = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate    = cur_rates[i]
    end
    # constant jump rates
    idx  += 1
    rates = p.rates
    @inbounds for i in 1:length(p.rates)
      new_rate       = rates[i](u, params, t)
      cur_rates[idx] = new_rate + prev_rate
      prev_rate      = cur_rates[idx]
      idx           += 1
    end

    @inbounds sum_rate = cur_rates[end]
    p.sum_rate = sum_rate
    nothing
end

function calculate_sum_rate(p::DelayDirectJumpAggregation, u,  params, t)
    prev_rate = zero(t)
    new_rate  = zero(t)
    cur_rates = zeros(length(p.cur_rates))
    majumps = p.ma_jumps
    idx       = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate     = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate    = cur_rates[i]
    end
    # constant jump rates
    idx  += 1
    rates = p.rates
    @inbounds for i in 1:length(p.rates)
      new_rate       = rates[i](u, params, t)
      cur_rates[idx] = new_rate + prev_rate
      prev_rate      = cur_rates[idx]
      idx           += 1
    end
    cur_rates[end]
end
"""
    function create_Tstruct(de_chan::Vector{Vector{T}})

calculate `Tstruct` according to the de_chan

# Arguments
- `de_chan::Vector{Vector{T}}`: where T is the tpye of `t`
"""
function create_Tstruct(de_chan::Vector{Vector{T}}) where {T}
    N = sum(length.(de_chan))
    Tstruct1 = Vector{T}(undef,N)
    Tstruct2 = Vector{Int64}(undef,N)
    k = 1
    @inbounds while k <= N   
        for i in eachindex(de_chan)
            for j in eachindex(de_chan[i])
                Tstruct1[k] = de_chan[i][j]
                Tstruct2[k] = i
                k += 1
            end
        end
    end
    vecorder = sortperm(Tstruct1)
    Tstruct1 = Tstruct1[vecorder]
    Tstruct2 = Tstruct2[vecorder]
    Tstruct1, Tstruct2
end
"""
Update state according up the next_jump;
"""
@inline function update_state_delay_Direct!(p::DelayDirectJumpAggregation, integrator, u, t)
    @unpack ma_jumps, next_jump, time_to_next_jump = p
    @unpack delay_interrupt, delay_trigger, delay_trigger_set, delay_interrupt_set = integrator.delayjumpsets

    ttnj = time_to_next_jump
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
    # shift delay channel !
    shift_delay_channel!(integrator.de_chan,ttnj) 
    update_delay_channel!(integrator.de_chan)
    if next_jump in delay_interrupt_set
        # delay_chan is changed according to affect_chan!
        delay_interrupt[next_jump](integrator.de_chan, p.rng) #affect_chan! decides how to modify the molecules in the delay channels
    elseif next_jump in delay_trigger_set
        # delay_trigger[next_jump] will affect de_channel
        delay_trigger[next_jump](integrator.de_chan, p.rng) # change to affect!
    end
    # save jump that was just executed
    p.prev_jump = next_jump
    nothing
end

