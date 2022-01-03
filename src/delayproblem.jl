
#BEGIN DelayJump
mutable struct DelayJumpProblem{iip,P,A,C,J<:Union{Nothing,AbstractJumpAggregator},J2,J3,J4,J5,deType} <: DiffEqBase.AbstractJumpProblem{P,J}
    prob::P
    aggregator::A
    discrete_jump_aggregation::J
    jump_callback::C
    variable_jumps::J2
    regular_jump::J3
    massaction_jump::J4
    delayjumpsets::J5
    de_chan0::deType
end

function DelayJumpProblem(jprob::JumpProblem, delayjumpsets::DelayJumpSet,de_chan0)
    @unpack prob, aggregator, discrete_jump_aggregation, jump_callback, variable_jumps, regular_jump, massaction_jump = jprob
    if !(aggregator<:AbstractDSSAJumpAggregator)
        error("To solve DelayJumpProblem, one has to use one of the delay aggregators.")
    end
    DelayJumpProblem(prob, aggregator, discrete_jump_aggregation, jump_callback, variable_jumps, regular_jump, massaction_jump, delayjumpsets, de_chan0)
end

function DelayJumpProblem(p::P,a::A,dj::J,jc::C,vj::J2,rj::J3,mj::J4,djs::J5,de_chan0::deType) where {P,A,J,C,J2,J3,J4,J5,deType}
    iip = isinplace_jump(p,rj)
    DelayJumpProblem{iip,P,A,C,J,J2,J3,J4,J5,deType}(p,a,dj,jc,vj,rj,mj,djs,de_chan0)
end


"""
    function DelayJumpProblem(prob::DiscreteProblem, aggregator::AbstractDelayAggregatorAlgorithm, jumps::JumpSet, delayjumpsets::DelayJumpSet, de_chan0)
# Fields
- `prob::DiscreteProblem`
- `aggregator::AbstractDelayAggregatorAlgorithm`
- `jumps::JumpSet`
- `delayjumpsets::DelayJumpSet`
- `de_chan0::Vector{Vector{T}}` 
TODO
"""
function DelayJumpProblem(prob, aggregator::AbstractDelayAggregatorAlgorithm, jumps::JumpSet, delayjumpsets::DelayJumpSet, de_chan0;
                     save_positions = typeof(prob) <: DiffEqBase.AbstractDiscreteProblem ? (false,true) : (true,true),
                     rng = Xorshifts.Xoroshiro128Star(rand(UInt64)), scale_rates = true, useiszero = true, spatial_system=nothing, hopping_constants=nothing, kwargs...)

  # initialize the MassActionJump rate constants with the user parameters
  if using_params(jumps.massaction_jump) 
    rates = jumps.massaction_jump.param_mapper(prob.p)
    maj = MassActionJump(rates, jumps.massaction_jump.reactant_stoch, jumps.massaction_jump.net_stoch, 
                         jumps.massaction_jump.param_mapper; scale_rates=scale_rates, useiszero=useiszero, 
                         nocopy=true)
  else
    maj = jumps.massaction_jump
  end

  ## Spatial jumps handling
  if spatial_system !== nothing && hopping_constants !== nothing && !is_spatial(aggregator) # check if need to flatten
    prob, maj = flatten(maj, prob, spatial_system, hopping_constants; kwargs...)
  end
  ## Constant Rate Handling
  t,end_time,u = prob.tspan[1],prob.tspan[2],prob.u0
  if (typeof(jumps.constant_jumps) <: Tuple{}) && (maj === nothing) && !is_spatial(aggregator) # check if there are no jumps
    disc = nothing
    constant_jump_callback = CallbackSet()
  else
    disc = aggregate(aggregator,u,prob.p,t,end_time,jumps.constant_jumps,maj,save_positions,rng; spatial_system = spatial_system, hopping_constants = hopping_constants, kwargs...)
    constant_jump_callback = DiscreteCallback(disc)
  end

  iip = isinplace_jump(prob, jumps.regular_jump)

  ## Variable Rate Handling
  if typeof(jumps.variable_jumps) <: Tuple{}
    new_prob = prob
    variable_jump_callback = CallbackSet()
  else
    new_prob = extend_problem(prob,jumps)
    variable_jump_callback = build_variable_callback(CallbackSet(),0,jumps.variable_jumps...)
  end
  callbacks = CallbackSet(constant_jump_callback,variable_jump_callback)

  DelayJumpProblem{iip,typeof(new_prob),typeof(aggregator),typeof(callbacks),
              typeof(disc),typeof(jumps.variable_jumps),
              typeof(jumps.regular_jump),typeof(maj),typeof(delayjumpsets),typeof(de_chan0)}(
                        new_prob,aggregator,disc,
                        callbacks,
                        jumps.variable_jumps,
                        jumps.regular_jump, maj, delayjumpsets, de_chan0)
end

"""
    function DelayJumpProblem(js::JumpSystem, prob, aggregator, delayjumpset, de_chan0; kwargs...)
# Fields
- `js::JumpSystem`: ...
- `prob::DiscreteProblem`: ...
- `aggregator::AbstractDelayAggregatorAlgorithm`
- `delayjumpsets::DelayJumpSet`
- `de_chan0::Vector{Vector{T}}` 
TODO
"""
function DelayJumpProblem(js::JumpSystem, prob, aggregator, delayjumpset, de_chan0; kwargs...)
    statetoid = Dict(value(state) => i for (i,state) in enumerate(states(js)))
    eqs       = equations(js)
    invttype  = prob.tspan[1] === nothing ? Float64 : typeof(1 / prob.tspan[2])

    # handling parameter substition and empty param vecs
    p = (prob.p isa DiffEqBase.NullParameters || prob.p === nothing) ? Num[] : prob.p

    majpmapper = JumpSysMajParamMapper(js, p; jseqs=eqs, rateconsttype=invttype)
    majs = isempty(eqs.x[1]) ? nothing : assemble_maj(eqs.x[1], statetoid, majpmapper)
    crjs = ConstantRateJump[assemble_crj(js, j, statetoid) for j in eqs.x[2]]
    vrjs = VariableRateJump[assemble_vrj(js, j, statetoid) for j in eqs.x[3]]
    ((prob isa DiscreteProblem) && !isempty(vrjs)) && error("Use continuous problems such as an ODEProblem or a SDEProblem with VariableRateJumps")
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, majs)

    if needs_vartojumps_map(aggregator) || needs_depgraph(aggregator)
        jdeps = asgraph(js)
        vdeps = variable_dependencies(js)
        vtoj = jdeps.badjlist
        jtov = vdeps.badjlist
        jtoj = needs_depgraph(aggregator) ? eqeq_dependencies(jdeps, vdeps).fadjlist : nothing
    else
        vtoj = nothing; jtov = nothing; jtoj = nothing
    end
    DelayJumpProblem(prob, aggregator, jset, delayjumpset, de_chan0; dep_graph=jtoj, vartojumps_map=vtoj, jumptovars_map=jtov, scale_rates=false, nocopy=true, kwargs...)
end


"""
for remaking
"""
function DiffEqBase.remake(thing::DelayJumpProblem; kwargs...)

  errmesg = """
  DelayJumpProblems can currently only be remade with new u0, de_chan0, p, tspan, delayjumpsets fields, prob fields. 
  """
  !issubset(keys(kwargs),((:u0,:de_chan0,:p,:tspan,:prob)...,propertynames(thing.delayjumpsets)...)) && error(errmesg)

  if :prob ∉ keys(kwargs)
    dprob = DiffEqBase.remake(thing.prob; kwargs...)
    # if the parameters were changed we must remake the MassActionJump too
    if (:p ∈ keys(kwargs)) && DiffEqJump.using_params(thing.massaction_jump)
        DiffEqJump.update_parameters!(thing.massaction_jump, dprob.p; kwargs...)
    end      
  else
    any(k -> k in keys(kwargs), (:u0,:p,:tspan)) && error("If remaking a DelayJumpProblem you can not pass both prob and any of u0, p, or tspan.")
    dprob = kwargs[:prob]

    # we can't know if p was changed, so we must remake the MassActionJump
    if DiffEqJump.using_params(thing.massaction_jump)
      DiffEqJump.update_parameters!(thing.massaction_jump, dprob.p; kwargs...)
    end 
  end
  if any(k -> k in keys(kwargs), propertynames(thing.delayjumpsets)) 
      delayjumpsets = update_delayjumpsets(thing.delayjumpsets; kwargs...)
  else
      delayjumpsets = thing.delayjumpsets
  end
  de_chan0 = :de_chan0 ∈ keys(kwargs) ? kwargs[:de_chan0] : thing.de_chan0

  DelayJumpProblem(dprob, thing.aggregator, thing.discrete_jump_aggregation, thing.jump_callback,
     thing.variable_jumps, thing.regular_jump, thing.massaction_jump, delayjumpsets, de_chan0)
end

function update_delayjumpsets(delayjumpsets::DelayJumpSet; kwargs...)
    delayjumpsets_ = deepcopy(delayjumpsets)
    for (key, value) in kwargs
        if toexpr(key) in [:delay_trigger,:delay_interrupt, :delay_complete] 
            setproperty!(delayjumpsets_, toexpr(key), value)
        end
    end    
    setproperty!(delayjumpsets_, :delay_interrupt_set, collect(keys(delayjumpsets_.delay_interrupt)))
    setproperty!(delayjumpsets_, :delay_trigger_set, collect(keys(delayjumpsets_.delay_trigger)))
    return delayjumpsets_
end

Base.summary(io::IO, prob::DelayJumpProblem) = string(DiffEqBase.parameterless_type(prob)," with problem ",DiffEqBase.parameterless_type(prob.prob)," and aggregator ",typeof(prob.aggregator))
function Base.show(io::IO, mime::MIME"text/plain", A::DelayJumpProblem)
  println(io,summary(A))
  println(io,"Number of constant rate jumps: ",A.discrete_jump_aggregation === nothing ? 0 : num_constant_rate_jumps(A.discrete_jump_aggregation))
  println(io,"Number of variable rate jumps: ",length(A.variable_jumps))
  if A.regular_jump !== nothing
    println(io,"Have a regular jump")
  end
  if (A.massaction_jump !== nothing) && (get_num_majumps(A.massaction_jump) > 0)
    println(io,"Have a mass action jump")
  end
  if A.delayjumpsets !== nothing
    println(io,"Number of delay trigger reactions: ",length(A.delayjumpsets.delay_trigger))
    println(io,"Number of delay interrupt reactions: ",length(A.delayjumpsets.delay_interrupt))
  end
end
#END DelayJump