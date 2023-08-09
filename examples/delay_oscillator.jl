using DelaySSAToolkit
using Catalyst
# using Plots
################################################################################
# the following example illustrates how to deal with a oscillatory system
# d: X-> 0; k1/(1+Y^2): 0-> X; [trigger X-> Y after τ time;] k2*Y/(1+Y):  Y-> 0;  

rn = @reaction_network begin
    1 / (1 + Y^2), 0 --> X
    1 / (1 + Y), Y --> 0
end

jumpsys = convert(JumpSystem, rn; combinatoric_ratelaws=false)

u0 = [0, 0]
de_chan0 = [[]]
tf = 400.0
tspan = (0, tf)
τ = 20.0
dprob = DiscreteProblem(jumpsys, u0, tspan)
# jumpsys.dep_graph

delay_trigger_affect! = function (integrator, rng)
    return append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1 => delay_trigger_affect!)
delay_complete = Dict(1 => [2 => 1, 1 => -1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
# alg = DelayRejection()
# alg = DelayDirect()
# alg = DelayMNRM()
alg = DelayDirectCR()
djprob = DelayJumpProblem(jumpsys, dprob, alg, delayjumpset, de_chan0)
sol = solve(djprob, SSAStepper(); seed=12345)
using Plots
plot(
    sol.t,
    [sol.u[i][1] for i in eachindex(sol.u)];
    alpha=0.3,
    label="X",
    color="green",
    linewidth=2,
    legend=:top,
    ylabel="# of individuals",
    xlabel="Time",
    fmt=:svg,
)

sol = solve(djprob, SSAStepper(); seed=1234)
plot!(
    sol.t,
    [sol.u[i][1] for i in eachindex(sol.u)];
    label="X",
    linewidth=2,
    line=:dash,
    color="green",
)

Sample_size = Int(1e4)
ens_prob = EnsembleProblem(djprob)
ens = solve(ens_prob, SSAStepper(), EnsembleThreads(); trajectories=Sample_size, saveat=0.1)
using StatsBase
tsm(t) = timepoint_mean(ens, t)
mean_X(t) = tsm(t)[1]
mean_Y(t) = tsm(t)[2]

timestamps = 0:1.0:tf
plot(
    timestamps,
    mean_X.(timestamps);
    linewidth=3,
    line=:dash,
    label="X",
    xlabel="time",
    ylabel="Mean Value",
)
plot!(timestamps, mean_Y.(timestamps); linewidth=3, line=:dash, legend=:topright, label="Y")

# plot(mean_X.(timestamps),mean_Y.(timestamps))
# might be very slow, can try with small end point
# @gif for i in 1:200
#     plot(mean_X,mean_Y,0,i, linewidth=range(0, 5, length = 200),seriesalpha=range(0, 1, length = 200),xlim=(0,21),ylim=(0,7),label=false,xlabel="X",ylabel="Y")
# end
