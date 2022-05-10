using Catalyst, DelaySSAToolkit
# 0 -->A1, => A2, ..., => AN => 0 with delay 1 for each cascade delay, N being the length of the delay reaction chain.

rn = @reaction_network begin
    5, 0 --> A1
end

delay_trigger_affect! = []
chain_len = 10
delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], 1.)
end


u0 = zeros(chain_len)
de_chan0 = []
for _ in 1:chain_len
    push!(de_chan0, [])
end
tspan = (0.,50.)
delay_complete_affect! = []
for i in 1:chain_len-1
    push!(delay_complete_affect!, function (integrator, rng)
        integrator.u[i] -= 1 # A_prev minus 1
        integrator.u[i+1] += 1 # A plus 1
        append!(integrator.de_chan[i+1], 1.) # add to the delay channel
    end
    )
end
push!(delay_complete_affect!, function (integrator, rng)
    integrator.u[chain_len] -= 1 # A_prev minus 1
end)

delay_trigger = Dict(Pair(1, delay_trigger_affect!))
delay_complete = Dict(i=>delay_complete_affect![i] for i in 1:chain_len)
delay_interrupt = Dict()


delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws = false)
dprob = DiscreteProblem(jumpsys, u0, tspan)


using Test
algos = [DelayDirect(), DelayRejection(), DelayMNRM(),  DelayDirectCR()]
@testset for algo in algos
    djprob = DelayJumpProblem(jumpsys, dprob, algo, delayjumpset, de_chan0, save_positions=(false,false), save_delay_channel= true)
    sol = solve(djprob, SSAStepper(), saveat = 1.)
    for i in 1:length(sol.t)-1, j in 1:chain_len-1 
        err = sum(abs2, sol.channel[i][j] .- sol.channel[i+1][j+1])
        @test err < 1e-20
    end        
end
