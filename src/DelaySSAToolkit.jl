module DelaySSAToolkit

using Reexport
@reexport using DiffEqBase
using UnPack, Random, LinearAlgebra, DiffEqBase, StaticArrays, DocStringExtensions
using FunctionWrappers, RandomNumbers
using Base.FastMath: add_fast, sub_fast, exp_fast, log_fast, div_fast

import DiffEqBase: DiscreteCallback, init, solve, solve!, initialize!
# Types and Structs
using JumpProcesses
import JumpProcesses: AbstractJumpAggregator, AbstractJump, JumpProblem, ConstantRateJump, VariableRateJump, RegularJump, MassActionJump, JumpSet, SSAStepper
# Dependency graph functions
import JumpProcesses: make_dependency_graph, add_self_dependencies!, var_to_jumps_map
# other functions
import JumpProcesses: using_params, get_jump_info_fwrappers, isinplace_jump, extend_problem, build_variable_callback, get_num_majumps

# using Catalyst 
using ModelingToolkit
import ModelingToolkit: JumpSysMajParamMapper, assemble_maj, assemble_vrj, assemble_crj, eqeq_dependencies, variable_dependencies, JumpSystem, states, equations, asgraph, value

using DataStructures
import DataStructures: update!


using StaticArrays, Base.Threads
import Base.Threads
@static if VERSION < v"1.3"
  seed_multiplier() = Threads.threadid()
else
  seed_multiplier() = 1
end

include("delayaggregator/aggregators.jl")
include("delayproblem.jl")
include("delayaggregator/delayssajump.jl")
include("delayaggregator/delayrejection.jl")
include("delayaggregator/delaymnrm.jl")
include("delayaggregator/delaydirect.jl")
include("delayaggregator/delaydirectCR.jl")
include("delayaggregator/delaycoevolve.jl")
include("delaySSA_stepper.jl")
include("utils.jl")

# export DSSAIntegrator
export DelayJumpProblem, DelayJumpSet, SSAStepper, MassActionJump, ConstantRateJump, JumpProblem, JumpSystem, JumpSet
export DelayRejection, DelayMNRM, DelayDirect, DelayDirectCR, DelayCoevolve
export solve, remake

end
