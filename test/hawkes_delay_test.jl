# this test is a small modification of https://github.com/SciML/JumpProcesses.jl/blob/master/test/hawkes_test.jl
using DelaySSAToolkit, JumpProcesses, OrdinaryDiffEq, Statistics
using Test
using Random
rng = MersenneTwister(12345)
delay = 20.0
function reset_history!(h; start_time=nothing)
    @inbounds for i in 1:length(h)
        h[i] = eltype(h)[]
    end
    return nothing
end

function empirical_rate(sol, agg)
    if typeof(agg) <: DelayCoevolve
        return (sol(sol.t[end]) - sol(sol.t[1] + delay)) / (sol.t[end] - sol.t[1] - delay)
    else
        return (sol(sol.t[end]) - sol(sol.t[1])) / (sol.t[end] - sol.t[1])
    end
end

function hawkes_rate(i::Int, g, h)
    function rate(u, p, t)
        λ, α, β = p
        x = zero(typeof(t))
        for j in g[i]
            for _t in reverse(h[j])
                λij = α * exp(-β * (t - _t))
                if λij ≈ 0
                    break
                end
                x += λij
            end
        end
        return λ + x
    end
    return rate
end

function hawkes_jump(i::Int, g, h, agg; uselrate=true)
    rate = hawkes_rate(i, g, h)
    urate = rate
    if uselrate
        lrate(u, p, t) = p[1]
        rateinterval =
            (u, p, t) -> begin
                _lrate = lrate(u, p, t)
                _urate = urate(u, p, t)
                return _urate == _lrate ? typemax(t) : 1 / (2 * _urate)
            end
    else
        lrate = nothing
        rateinterval = (u, p, t) -> begin
            _urate = urate(u, p, t)
            return 1 / (2 * _urate)
        end
    end
    if typeof(agg) <: DelayCoevolve
        affect! = (integrator) -> begin
            push!(h[i], integrator.t)
            integrator.u[i] += 0
        end
    else
        affect! = (integrator) -> begin
            push!(h[i], integrator.t)
            integrator.u[i] += 1
        end
    end
    return VariableRateJump(rate, affect!; lrate, urate, rateinterval)
end

function hawkes_jump(u, g, h, agg; uselrate=true)
    return [hawkes_jump(i, g, h, agg; uselrate) for i in 1:length(u)]
end

function hawkes_problem(
    p,
    agg::DelayCoevolve;
    u=[0.0],
    tspan=(0.0, 50.0),
    save_positions=(false, true),
    g=[[1]],
    h=[[]],
    uselrate=true,
)
    dprob = DiscreteProblem(u, tspan, p)
    jumps = JumpSet(hawkes_jump(u, g, h, agg; uselrate)...)
    de_chan0 = [[]]
    delay_trigger = Dict(1 => [1 => delay]) # add a delay of 1.0 to the first jump
    delay_complete = Dict(1 => [1 => 1]) # complete the delay will duplicate 1 product
    delay_interrupt = Dict()
    delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
    jprob = DelayJumpProblem(
        dprob, agg, jumps, delayjumpset, de_chan0; dep_graph=g, save_positions, rng
    )
    return jprob
end

function f!(du, u, p, t)
    du .= 0
    return nothing
end

function hawkes_problem(
    p,
    agg;
    u=[0.0],
    tspan=(0.0, 50.0),
    save_positions=(false, true),
    g=[[1]],
    h=[[]],
    kwargs...,
)
    oprob = ODEProblem(f!, u, tspan, p)
    jumps = hawkes_jump(u, g, h, agg)
    jprob = JumpProblem(oprob, agg, jumps...; save_positions, rng)
    return jprob
end

function expected_stats_hawkes_problem(p, tspan, agg)
    if typeof(agg) <: DelayCoevolve
        T = tspan[end] - tspan[1] + delay
        # stepper = SSAStepper()
    else
        T = tspan[end] - tspan[1]
    end
    λ, α, β = p
    γ = β - α
    κ = β / γ
    Eλ = λ * κ
    # Equation 21
    # J. Da Fonseca and R. Zaatour,
    # “Hawkes Process: Fast Calibration, Application to Trade Clustering and Diffusive Limit.”
    # Rochester, NY, Aug. 04, 2013. doi: 10.2139/ssrn.2294112.
    Varλ = (Eλ * (T * κ^2 + (1 - κ^2) * (1 - exp(-T * γ)) / γ)) / (T^2)
    return Eλ, Varλ
end

u0 = [0.0]
p = (0.5, 0.5, 2.0)
tspan = (0.0, 250.0)
g = [[1]]
h = [Float64[]]

aggs = (Direct(), DelayCoevolve(), DelayCoevolve())
uselrate = zeros(Bool, length(aggs))
uselrate[3] = true
Nsims = Int(5e2)

for (i, agg) in enumerate(aggs)
    @info "Testing $(typeof(agg))"
    jump_prob = hawkes_problem(p, agg; u=u0, tspan, g, h, uselrate=uselrate[i])
    if typeof(agg) <: DelayCoevolve
        stepper = SSAStepper()
    else
        stepper = Tsit5()
    end
    sols = Vector{ODESolution}(undef, Nsims)
    for n in 1:Nsims
        reset_history!(h)
        sols[n] = solve(jump_prob, stepper)
    end
    if typeof(agg) <: DelayCoevolve
        λs = permutedims(mapreduce((sol) -> empirical_rate(sol, agg), hcat, sols))
    else
        cols = length(sols[1].u[1].u)
        λs = permutedims(mapreduce((sol) -> empirical_rate(sol, agg), hcat, sols))[
            :, 1:cols
        ]
    end
    Eλ, Varλ = expected_stats_hawkes_problem(p, tspan, agg)
    @test isapprox(mean(λs), Eλ; atol=0.01)
    @test isapprox(var(λs), Varλ; atol=0.001)
end
