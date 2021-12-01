using DelaySSAToolkit, Catalyst
using DiffEqJump

#1st route Reaction network
rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r
# Define Jumpsystem
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)


# Define DelayJumpProblem = DiscreteProblem + JumpSystem
u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2]
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
τ = 20.
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
# algo = DelayDirect()
algo = DelayRejection()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
jprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0, save_positions=(true,true))


# Solve
sol = solve(jprob, SSAStepper(), seed = 1234)
using Plots
fig = plot(sol, label = ["S" "I" "E" "R"], linewidth = 3, legend = :top, ylabel = "# of individuals", xlabel = "Time", fmt=:svg)
savefig(fig,"docs/src/assets/seir.svg")


# 2nd route Define DiscreteProblem + JumpSet
# S I E R 1,2,3,4
# ρ: S+I -> E+I
# r: I -> R
ρ, r = [1e-4, 1e-2]
rate1 = [ρ, r]
reactant_stoich = [[1=>1,2=>1],[2=>1]]
net_stoich = [[1=>-1,3=>1],[2=>-1,4=>1]]
mass_jump = MassActionJump(rate1, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,[mass_jump])


# Define DelayJumpProblem =  DiscreteProblem + JumpSet
u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
dprob = DiscreteProblem(u0, tspan)
τ = 20.
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
# algo = DelayDirect()
algo = DelayRejection()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
djprob = DelayJumpProblem(dprob, algo, jumpset, delayjumpset, de_chan0, save_positions=(true,true))

sol = solve(djprob, SSAStepper(), seed = 1234)
using Plots
fig = plot(sol, label = ["S" "I" "E" "R"], linewidth = 3, legend = :top, ylabel = "# of individuals", xlabel = "Time", fmt=:svg)


using PkgTemplates
tpl = Template(;
    user = "palmtree2013",
    authors = "Xiaoming Fu",
    julia=v"1.6",
    dir="~/.julia/dev",
    plugins=[
        Git(; manifest=true, ssh=false),
        GitHubActions(; x86=true),
        Codecov(),
        Documenter{GitHubActions}(),
        Develop(),
    ],
)

tpl("DelaySSAdocs.jl")