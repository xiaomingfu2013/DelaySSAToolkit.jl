using DelaySSAToolkit
using Test, SafeTestsets

@time begin
    @time @safetestset "Algorithm accuracy" begin
        include("bursty_model.jl")
    end 
    @time @safetestset "dep_gr_delay test" begin
        include("dep_gr_delay.jl")
    end 
    @time @safetestset "save delay channel test" begin
        include("save_delay_channel.jl")
    end 
    @time @safetestset "cascade of delay reaction test" begin
        include("cascade_of_delay_reaction.jl")
    end
    @time @safetestset "delay problem test" begin
        include("delay_problem_test.jl")
    end  
    @time @safetestset "remake problem test" begin
        include("remake_test.jl")
    end
    @time @safetestset "low level interface test" begin
        include("low_level_interface.jl")
    end
end
