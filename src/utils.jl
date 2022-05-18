



Base.@propagate_inbounds Base.getindex(A::DSSASolution, i::Int) = [A.u[i], A.channel[i]]
Base.@propagate_inbounds Base.getindex(A::DSSASolution, i::Int, ::Colon) = [A.u[j][i] for j in 1:length(A.t)]

function (A::DSSASolution)(s::Symbol,i::Int)
    if s == :channel
        @assert i <= length(A.channel[1])
        return [A.channel[j][i] for j in 1:length(A.t)]
    elseif s ==:u
        @assert i <= length(A.u[1])
        return [A.u[j][i] for j in 1:length(A.t)]
    end
end 

function (A::DSSASolution)(tval)
    @unpack odesol = A
    odesol(tval)
end 



(integrator::DSSAIntegrator)(t) = copy(integrator.u)
(integrator::DSSAIntegrator)(out,t) = (out .= integrator.u)


function Base.show(io::IO, m::MIME"text/plain", A::DSSASolution)
    println(io,string("retcode: ",A.odesol.retcode))
    println(io,string("Interpolation: "),DiffEqBase.interp_summary(A.odesol.interp))
    print(io,"t: ")
    show(io,m,A.t)
    println(io)
    print(io,"u: ")
    show(io,m,A.u)
    println(io)
    print(io,"channel: ")
    show(io,m,A.channel)
    println(io,"\n===\nUse sol.u to check the state variable and sol.channel to check the delay channel solution.\n===")
end



# """
#     function get_reaction_idx(rn::ReactionSystem)

# Get the rearranged order for a reaction system, returns a matrix `m` with first column as the rearranged reactions, second column is the order in the orign `ReactionSystem`
# If you know the old idx for a certain reaction in `ReactionSystem`, you can use 
# `m[:,2][idx]` to get the new idx. 
# """

# function get_reaction_idx(rn::Catalyst.ReactionSystem)
#     massactionjump_order = Int64[]
#     constantjump_order = Int64[]
#     rxvars = []
#     for (i, rx) in enumerate(ModelingToolkit.get_eqs(rn))
#         empty!(rxvars)
#         (Catalyst.reactionrates(rn) isa ModelingToolkit.Symbolic) && ModelingToolkit.get_variables!(rxvars, Catalyst.reactionrates(rn))
#         @inbounds for j = 1:length(rxvars)
#             if isequal(rxvars[j], ModelingToolkit.get_iv(rn))
#                 error("Does not support VariableRateJump")
#             end
#         end
#         if Catalyst.ismassaction(rx, rn)
#             push!(massactionjump_order, i)
#         else
#             push!(constantjump_order, i)
#         end
#     end
#     order = vcat(massactionjump_order, constantjump_order)
#     # println("The rearranged order :")
#     vcat([[ModelingToolkit.get_eqs(rn)[order[i]] order[i]] for i in eachindex(order)]...) 
# end
# export get_reaction_idx