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
    # @time @safetestset "Index test" begin
    #     include("check_index_reactions.jl")
    # end
end
