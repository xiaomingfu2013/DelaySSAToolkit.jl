export DSSAIntegrator


mutable struct DSSAIntegrator{F,uType,tType,P,S,CB,SA,OPT,TS,chanType,chanS} <: DiffEqBase.DEIntegrator{SSAStepper,Nothing,uType,tType}
    f::F
    u::uType
    t::tType
    tprev::tType
    p::P
    sol::S
    i::Int
    tstop::tType
    cb::CB
    saveat::SA
    save_everystep::Bool
    save_end::Bool
    cur_saveat::Int
    opts::OPT
    tstops::TS
    tstops_idx::Int
    u_modified::Bool
    keep_stepping::Bool  # false if should terminate a simulation
    de_chan::chanType
    chan_sol::chanS
    delayjumpsets::DelayJumpSet
    save_delay_channel::Bool
end


mutable struct DSSASolution{uType,tType,chanS}
    u::uType
    t::tType
    channel::chanS
    odesol::DiffEqBase.ODESolution
end
function DSSASolution(ode_sol::DiffEqBase.ODESolution, channel)
    @unpack u, t = ode_sol
    DSSASolution{typeof(u),typeof(t),typeof(channel)}(u, t, channel, ode_sol)
end




function DiffEqBase.u_modified!(integrator::DSSAIntegrator, bool::Bool)
    integrator.u_modified = bool
end

function DiffEqBase.__solve(djump_prob::DelayJumpProblem,
    alg::SSAStepper;
    kwargs...)
    integrator = init(djump_prob, alg; kwargs...)
    solve!(integrator)
    if integrator.save_delay_channel
        DSSASolution(integrator.sol, integrator.chan_sol)
    else
        integrator.sol
    end
end

## Initiate the Problem

function DiffEqBase.__init(djump_prob::DelayJumpProblem,
    alg::SSAStepper;
    save_start=true,
    save_end=true,
    seed=nothing,
    alias_jump=Threads.threadid() == 1,
    saveat=nothing,
    callback=nothing,
    tstops=eltype(djump_prob.prob.tspan)[],
    numsteps_hint=100)
    if !(djump_prob.prob isa DiscreteProblem)
        error("SSAStepper only supports DiscreteProblems.")
    end
    @assert isempty(djump_prob.jump_callback.continuous_callbacks) # check if continuous_callbacks is empty
    ## Set up the cb which encodes the aggregator information
    if alias_jump
        cb = djump_prob.jump_callback.discrete_callbacks[end]
        if seed !== nothing
            Random.seed!(cb.condition.rng, seed)
        end
    else
        cb = deepcopy(djump_prob.jump_callback.discrete_callbacks[end])
        if seed === nothing
            Random.seed!(cb.condition.rng, seed_multiplier() * rand(UInt64))
        else
            Random.seed!(cb.condition.rng, seed)
        end
    end

    opts = (callback=CallbackSet(callback),)
    prob = djump_prob.prob
    de_chan0 = convert(Vector{Vector{typeof(prob.tspan[1])}}, djump_prob.de_chan0)
    if save_start
        t = [prob.tspan[1]]
        u = [copy(prob.u0)]
        chan_sol = [deepcopy(de_chan0)] # DelaySSA: build chan_sol
    else
        t = typeof(prob.tspan[1])[]
        u = typeof(prob.u0)[]
        chan_sol = typeof(de_chan0)[] # DelaySSA: build chan_sol
    end

    sol = DiffEqBase.build_solution(prob, alg, t, u, dense=false, calculate_error=false, destats=DiffEqBase.DEStats(0), interp=DiffEqBase.ConstantInterpolation(t, u))
    save_everystep = any(cb.save_positions)

    if typeof(saveat) <: Number
        _saveat = prob.tspan[1]:saveat:prob.tspan[2]
    else
        _saveat = saveat
    end

    if _saveat !== nothing && !isempty(_saveat) && _saveat[1] == prob.tspan[1]
        cur_saveat = 2
    else
        cur_saveat = 1
    end

    if _saveat !== nothing && !isempty(_saveat)
        sizehint!(u, length(_saveat) + 1)
        sizehint!(t, length(_saveat) + 1)
    elseif save_everystep
        sizehint!(u, numsteps_hint)
        sizehint!(t, numsteps_hint)
    else
        sizehint!(u, save_start + save_end)
        sizehint!(t, save_start + save_end)
    end

    integrator = DSSAIntegrator(prob.f, copy(prob.u0), prob.tspan[1], prob.tspan[1], prob.p, sol, 1, prob.tspan[1], cb, _saveat, save_everystep, save_end, cur_saveat, opts, tstops, 1, false, true, deepcopy(de_chan0), chan_sol, djump_prob.delayjumpsets, djump_prob.save_delay_channel)

    cb.initialize(cb, integrator.u, prob.tspan[1], integrator)
    DiffEqBase.initialize!(opts.callback, integrator.u, prob.tspan[1], integrator)
    integrator
end

function DiffEqBase.add_tstop!(integrator::DSSAIntegrator, tstop)
    if tstop > integrator.t
        future_tstops = @view integrator.tstops[integrator.tstops_idx:end]
        insert_index = integrator.tstops_idx + searchsortedfirst(future_tstops, tstop) - 1
        Base.insert!(integrator.tstops, insert_index, tstop)
    end
end


function DiffEqBase.solve!(integrator::DSSAIntegrator)

    end_time = integrator.sol.prob.tspan[2]
    while should_continue_solve(integrator) # It stops before adding a tstop over
        step!(integrator)
    end
    last_t = integrator.t
    integrator.t = end_time

    if typeof(integrator.cb.affect!) <: DelayDirectJumpAggregation
        saveat_function_direct_method_test!(integrator, last_t)
    else
        saveat_function!(integrator, last_t)
    end

    if integrator.save_end && integrator.sol.t[end] != end_time
        saveat_end_function!(integrator, last_t)
    end
    DiffEqBase.finalize!(integrator.opts.callback, integrator.u, integrator.t, integrator)
end


function saveat_end_function!(integrator, prev_t)
    # save last t
    end_time = integrator.sol.prob.tspan[2]
    push!(integrator.sol.t, end_time)
    t_final_gap = end_time - prev_t

    # save last u 
    if typeof(integrator.cb.affect!) <: DelayDirectJumpAggregation
        update_delay_at_tstop_test!(integrator.cb.affect!, integrator, integrator.p, prev_t, t_final_gap)
    end
    push!(integrator.sol.u, copy(integrator.u))

    # save last de_chan
    if integrator.save_delay_channel
        last_chan = deepcopy(integrator.de_chan)
        if !(typeof(integrator.cb.affect!) <: DelayDirectJumpAggregation)
            shift_delay_channel!(last_chan, t_final_gap)
            update_delay_channel!(last_chan)
        end
        push!(integrator.chan_sol, last_chan)
    end
end

"""
    function saveat_function!(integrator, prev_t)

Note that this function does not change u and de_chan, but only takes the snapshots of u and de_chan accordingly
"""
function saveat_function!(integrator, prev_t)
    if integrator.saveat !== nothing && !isempty(integrator.saveat)
        # Split to help prediction
        last_saved_t = prev_t
        while integrator.cur_saveat < length(integrator.saveat) && (integrator.saveat[integrator.cur_saveat] < integrator.t)
            tgap = integrator.saveat[integrator.cur_saveat] - last_saved_t
            push!(integrator.sol.t, integrator.saveat[integrator.cur_saveat])
            push!(integrator.sol.u, copy(integrator.u))
            if integrator.save_delay_channel
                prev_de_chan = deepcopy(integrator.de_chan)
                shift_delay_channel!(prev_de_chan, tgap)
                push!(integrator.chan_sol, prev_de_chan)
            end
            last_saved_t = integrator.saveat[integrator.cur_saveat]
            integrator.cur_saveat += 1
        end
    end
end

"""
    function saveat_function_direct_method_test!(integrator, prev_t)

Note that this function does not change u and de_chan for the next jump, but only takes the snapshots of u and de_chan accordingly. Because in next `execute_jumps`, shadow_integrator will overwrite integrator and when t = tend, the initial state will be recovered
"""
function saveat_function_direct_method_test!(integrator, prev_t)
    # Special to Direct method
    if integrator.saveat !== nothing && !isempty(integrator.saveat)
        last_saved_t = prev_t
        # save for last step and copied = false
        changed = false
        aggregator = integrator.cb.affect!
        if (integrator.cur_saveat < length(integrator.saveat)) && (integrator.saveat[integrator.cur_saveat] < integrator.t)
            shadow_integrator = aggregator.shadow_integrator
            shadow_integrator.u = copy(integrator.u)
            shadow_integrator.de_chan = deepcopy(integrator.de_chan)
            aggregator.copied = true
            changed = true
        end
        while integrator.cur_saveat < length(integrator.saveat) && (integrator.saveat[integrator.cur_saveat] < integrator.t)

            tgap = integrator.saveat[integrator.cur_saveat] - last_saved_t
            push!(integrator.sol.t, integrator.saveat[integrator.cur_saveat])

            update_delay_at_tstop_test!(aggregator, shadow_integrator, integrator.p, last_saved_t, tgap)
            push!(integrator.sol.u, copy(shadow_integrator.u))

            if integrator.save_delay_channel
                prev_de_chan = deepcopy(shadow_integrator.de_chan)
                push!(integrator.chan_sol, prev_de_chan)
            end

            last_saved_t = integrator.saveat[integrator.cur_saveat]
            integrator.cur_saveat += 1
        end
        if aggregator.copied && changed
            tend = integrator.t == integrator.sol.prob.tspan[2] ? integrator.sol.prob.tspan[2] : aggregator.next_jump_time
            last_tgap = tend - last_saved_t
            update_delay_at_tstop_test!(aggregator, shadow_integrator, integrator.p, last_saved_t, last_tgap)
        end
    end
end

# The Jump aggregators should not register the next jump through add_tstop! for SSAIntegrator
# such that we can achieve maximum performance
@inline function register_next_jump_time!(integrator::DSSAIntegrator, p::AbstractDSSAJumpAggregator, t)
    integrator.tstop = p.next_jump_time
    nothing
end


## Delay SSA modify
function DiffEqBase.step!(integrator::DSSAIntegrator)
    integrator.tprev = integrator.t
    next_jump_time = integrator.tstop > integrator.t ? integrator.tstop : typemax(integrator.tstop)

    doaffect = false
    if !isempty(integrator.tstops) &&
       integrator.tstops_idx <= length(integrator.tstops) &&
       integrator.tstops[integrator.tstops_idx] < next_jump_time

        integrator.t = integrator.tstops[integrator.tstops_idx]
        integrator.tstops_idx += 1
    else
        integrator.t = integrator.tstop
        doaffect = true # delay effect until after saveat
    end

    if typeof(integrator.cb.affect!) <: DelayDirectJumpAggregation
        saveat_function_direct_method_test!(integrator, integrator.tprev)
    else
        saveat_function!(integrator, integrator.tprev)
    end


    # FP error means the new time may equal the old if the next jump time is 
    # sufficiently small, hence we add this check to execute jumps until
    # this is no longer true.
    while integrator.t == integrator.tstop
        doaffect && integrator.cb.affect!(integrator)
    end

    if !(typeof(integrator.opts.callback.discrete_callbacks) <: Tuple{})
        discrete_modified, saved_in_cb = DiffEqBase.apply_discrete_callback!(integrator, integrator.opts.callback.discrete_callbacks...)
    else
        saved_in_cb = false
    end

    !saved_in_cb && savevalues!(integrator)

    nothing
end


function DiffEqBase.savevalues!(integrator::DSSAIntegrator, force=false)
    saved, savedexactly = false, false

    # No saveat in here since it would only use previous values,
    # so in the specific case of DSSAStepper it's already handled

    if integrator.save_everystep || force
        saved = true
        savedexactly = true
        push!(integrator.sol.t, integrator.t)
        push!(integrator.sol.u, copy(integrator.u))
        if integrator.save_delay_channel
            push!(integrator.chan_sol, deepcopy(integrator.de_chan))
        end
    end

    saved, savedexactly
end

function should_continue_solve(integrator::DSSAIntegrator)
    end_time = integrator.sol.prob.tspan[2]

    # we continue the solve if there is a tstop between now and end_time
    has_tstop = !isempty(integrator.tstops) &&
                integrator.tstops_idx <= length(integrator.tstops) &&
                integrator.tstops[integrator.tstops_idx] < end_time

    # we continue the solve if there will be a jump between now and end_time
    has_jump = integrator.t < integrator.tstop < end_time

    integrator.keep_stepping && (has_jump || has_tstop)
end


num_constant_rate_jumps(aggregator::AbstractDSSAJumpAggregator) = length(aggregator.rates)
