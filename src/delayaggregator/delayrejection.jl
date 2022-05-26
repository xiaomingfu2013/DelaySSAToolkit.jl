mutable struct DelayRejectionJumpAggregation{T,S,F1,F2,RNG} <: AbstractDSSAJumpAggregator
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
    next_delay::Vector{Int}
    num_next_delay::Vector{Int}
    time_to_next_jump::T
    dt_delay::T
end
#3
function DelayRejectionJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; kwargs...) where {T,S,F1,F2,RNG}
    nd = Int64[]
    nnd = Int64[]
    ttnj = zero(et)
    dt_delay = zero(et)
    #4
    DelayRejectionJumpAggregation{T,S,F1,F2,RNG}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, nd, nnd, ttnj, dt_delay)
end
#1 
function aggregate(aggregator::DelayRejection, u, p, t, end_time, constant_jumps,
    ma_jumps, save_positions, rng; kwargs...)


    # handle constant jumps using function wrappers
    rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)
    build_jump_aggregation(DelayRejectionJumpAggregation, u, p, t, end_time, ma_jumps, rates, affects!, save_positions, rng; kwargs...)
end



function initialize!(p::DelayRejectionJumpAggregation, integrator, u, params, t)
    find_next_delay_dt!(p, integrator)
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

@inline function execute_jumps!(p::DelayRejectionJumpAggregation, integrator, u, params, t)
    # update_state!(p, integrator, u)
    update_state_delay!(p, integrator, u, t)
    nothing
end

function generate_jumps!(p::DelayRejectionJumpAggregation, integrator, u, params, t)
    time_to_next_jump!(p, integrator, u, params, t)
    if isempty(p.next_delay)
        @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rand(p.rng) * p.sum_rate) # cur_rates is the cumsum of rates
    else
        p.next_jump = 0
    end
    nothing
end

@fastmath function time_to_next_jump!(p::DelayRejectionJumpAggregation, integrator, u, params, t)
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
    @inbounds for i in eachindex(p.rates)
        new_rate = rates[i](u, params, t)
        cur_rates[idx] = new_rate + prev_rate
        prev_rate = cur_rates[idx]
        idx += 1
    end

    @inbounds sum_rate = cur_rates[end]
    dt_reaction = randexp(p.rng) / sum_rate
    # consider delay dt 
    dt_delay_generation!(p, integrator)
    compare_delay!(p, integrator.de_chan, p.dt_delay, dt_reaction, t)
    p.sum_rate = sum_rate
    nothing
end


