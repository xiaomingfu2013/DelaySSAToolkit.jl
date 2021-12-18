using DelaySSAToolkit
using Test, SafeTestsets

@time begin
    @time @safetestset "Algorithm accuracy" begin
        include("bursty_model.jl")
    end 
end
