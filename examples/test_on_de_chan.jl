"""
find the minimal dt_delay in various delay channel
"""

function find_next_delay_dt(de_chan::Vector{Vector{T}}) where {T}
    val_vec = vcat(de_chan...)
    isempty(val_vec) ? Inf : minimum(val_vec)
end

function find_next_delay(de_chan::Vector{Vector{T}}) where {T}
    min_vec = Vector{T}(undef,length(de_chan))
    for i in eachindex(de_chan)
        min_vec[i] = isempty(de_chan[i]) ? Inf : minimum(de_chan[i])
    end
    findmin(min_vec)
end

# for vec
function find_next_delay_vec(de_chan::Vector{Vector{T}}, dt::T) where {T}
    position_vec = findall.(isequal(dt),de_chan)
    next_delay_vec =findall(!isempty,position_vec)
    num_next_delay_vec = length.(position_vec[next_delay_vec])
    return next_delay_vec, num_next_delay_vec
end

t = 0.
de_chan = [[0.09;0.09;rand(100)],[0.09,0.3],[0.1]]


de_chan = convert(Vector{Vector{typeof(t)}},de_chan)

# @benchmark vcat(de_chan...)

find_next_delay_dt(de_chan)
find_next_delay(de_chan)
find_next_delay_vec(de_chan, 0.09)

using BenchmarkTools
@benchmark find_next_delay_dt(de_chan)
@benchmark find_next_delay(de_chan)
@benchmark find_next_delay_vec(de_chan, 0.09)


@benchmark reduce(vcat,de_chan)
@benchmark collect(Iterators.flatten(de_chan))
# de_chan

# reduce(min,de_chan)