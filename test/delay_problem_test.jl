using Test, DelaySSAToolkit, Catalyst

ρ, r = [1e-4, 1e-2]
rates = [ρ, r]
reactant_stoich = [[1=>1,2=>1],[2=>1]]
net_stoich = [[1=>-1,3=>1],[2=>-1,4=>1]]
mass_jump = MassActionJump(rates, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,mass_jump)

u0 = [999,1,0,0] # S, I, E, R
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2] # parameters for ρ, r
τ = 20.
dprob = DiscreteProblem(u0,tspan,ps)

delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
de_chan0 = [[]]

djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions=(false, false), save_delay_channel = true)
sol = solve(djprob, SSAStepper(), seed = 1, saveat =10.)

io = IOBuffer()
show(io, "text/plain", djprob)
@test read(seekstart(io), String) == "\nNumber of constant rate jumps: 0\nNumber of variable rate jumps: 0\nHave a mass action jump\nNumber of delay trigger reactions: 1\nNumber of delay interrupt reactions: 0\n"

sol1 = solve(djprob, SSAStepper(), seed = 1, save_start = false)
io1 = IOBuffer()
show(io1, "text/plain", sol1)
@test read(seekstart(io1), String) == "retcode: Default\nInterpolation: Piecewise constant interpolation\nt: 1-element Vector{Float64}:\n 400.0\nu: 1-element Vector{Vector{Int64}}:\n [0, 129, 0, 871]\nchannel: 1-element Vector{Vector{Vector{Float64}}}:\n [[]]\n===\nUse sol.u to check the state variable and sol.channel to check the delay channel solution.\n===\n"


rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
djprob_ = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(false,false), save_delay_channel = true)

sol_ = solve(djprob_, SSAStepper(), seed =1, saveat = 10.)

@testset for i in eachindex(sol.t)
    @test sol[i] == sol_[i]
end

@testset for i in 1:3
    @test sol[i,:] == sol_[i,:]
    @test sol(:u, i) == sol_(:u, i)
end
@testset for t in sol.t
    @test sol(t) == sol_(t)
end
@test sol(:channel, 1) == sol_(:channel, 1)