module DelaySSAToolkit

using Reexport
@reexport using DiffEqBase
import DiffEqBase: DiscreteCallback, init, solve, solve!, plot_indices, initialize!
# Types and Structs
import DiffEqJump: AbstractAggregatorAlgorithm, AbstractJumpAggregator, AbstractJump, AbstractSSAJumpAggregator, JumpProblem, ConstantRateJump, VariableRateJump, RegularJump, MassActionJump, JumpSet, SSAStepper
# Dependency graph functions
import DiffEqJump: make_dependency_graph, add_self_dependencies!, var_to_jumps_map
# other functions
import DiffEqJump: using_params, get_jump_info_fwrappers, build_jump_aggregation, isinplace_jump, extend_problem, build_variable_callback, get_num_majumps, evalrxrate, executerx!, executerx

import Catalyst: JumpSystem, states, equations, assemble_jumps, asgraph,  variable_dependencies, eqeq_dependencies, value
import ModelingToolkit: JumpSysMajParamMapper, assemble_maj, assemble_vrj, assemble_crj
import DataStructures: update!

using DataStructures, UnPack, Random, LinearAlgebra, DiffEqBase, StaticArrays
using FunctionWrappers, ArrayInterface, TreeViews, RandomNumbers
import RecursiveArrayTools: recursivecopy!
using StaticArrays, Base.Threads
import Base.Threads
@static if VERSION < v"1.3"
  seed_multiplier() = Threads.threadid()
else
  seed_multiplier() = 1
end
include("delaysystem.jl")
include("delayaggregator/aggregators.jl")
include("delayaggregator/delayssajump.jl")
include("delayaggregator/delayrejection.jl")
include("delayaggregator/delaymnrm.jl")
include("delayaggregator/delaydirect.jl")

include("delayproblem.jl")

include("delaySSA_stepper.jl")
export DelayJumpProblem, DelaySSAsolve, DelayJumpSet
export DelayRejection, DelayMNRM, DelayDirect
export solve

end
