using Test, DelaySSAToolkit, Catalyst, StaticArrays

params = [0.1, 0.1, 10.0, 0.1, 10.0, 50.0]
σ_off, σ_on, ρ_on, d, τ, tf = params
u0_ =@SVector [0,1,0]
tspan = (0., 50.)


rates_ = [σ_off, σ_on, ρ_on]
dprob_ = DiscreteProblem(u0_, tspan, rates_)
delay_trigger = Dict(3 => [1 => τ])
delay_complete = Dict(1 => [3 => -1])
delay_affect! = function (integrator, rng)
    i = rand(rng, 1:length(integrator.de_chan[1]))
    deleteat!(integrator.de_chan[1], i)
end
delay_interrupt = Dict(4 => delay_affect!)

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)



rates_ = [σ_off, σ_on, ρ_on]
react_stoich_ = [[1 => 1], [2 => 1], [1 => 1]]
net_stoich_ = [[1 => -1, 2 => 1], [1 => 1, 2 => -1], [3 => 1]]
mass_action_jump_ = MassActionJump(rates_, react_stoich_, net_stoich_; scale_rates=false)

affect! = function (integrator)
    # integrator.u[3] -= 1
    integrator.u = setindex(integrator.u, integrator.u[3] - 1, 3)
end
rate2 =  (u,p,t) -> 0.1*u[3]
constant_rate_jump = ConstantRateJump(rate2, affect!)
jumpset_ = JumpSet((), (constant_rate_jump,), nothing, mass_action_jump_)



dep_gr_delay = Dict(1=>[4])
dep_gr = [[1,2,3],[1,2,3],[3,4],[4]]
de_chan0 = [[]]
algs = [DelayRejection(), DelayDirect(), DelayMNRM(), DelayDirectCR()]
alg = algs[4]


for alg in algs
    djprob__ = DelayJumpProblem(dprob_, alg, jumpset_, delayjumpset, de_chan0, save_positions=(true, true), save_delay_channel=true, dep_graph = dep_gr, dep_graph_delay = dep_gr_delay)
    @info "Testing method $(alg)"
    sol = solve(djprob__, SSAStepper(), seed =1)
    
    djprob_no_dep_gr_delay = DelayJumpProblem(dprob_, alg, jumpset_, delayjumpset, de_chan0, save_positions=(true, true), save_delay_channel=true, dep_graph = dep_gr)
    @info "Testing method $(alg)"
    sol_no_dep_gr_delay = solve(djprob_no_dep_gr_delay, SSAStepper(), seed = 1)
    
    for i in eachindex(sol.u)
        @test sol.u[i][3] == length(sol.channel[i][1])
    end

    for i in eachindex(sol.u)
        @test sol_no_dep_gr_delay.u[i][3] == length(sol_no_dep_gr_delay.channel[i][1])
    end
    for i in eachindex(sol.u)
        @test sol_no_dep_gr_delay.u[i] == sol.u[i]
        @test sol_no_dep_gr_delay.channel[i] == sol.channel[i]
    end
end

