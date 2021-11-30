export find_next_delay_vec
"""
    function isequal_in_vec(A::Vector{T}, x::T)

Find if there is a floating point number `x` in the vector `vec`.If `x` exists, return true; if not, return false.

# Examples
```julia-repl
julia> isequal_vec([1,2,3]], 1)
true
```
"""
@inline function isequal_in_vec(A::Vector{T}, x::T) where {T}
    count = 0
    for i in eachindex(A)
        @inbounds if A[i] == x
            count += 1
        end
    end
    if !isequal(0,count)
        return true
    else
        return false
    end
end

"""
    function find_num_in_vec(A::Vector{Vector{T}}, position_index::Vector{Int64}, x::T)

Find the number of values which in each vector elements equal to `x` according to the corresponding index position specified by the element in the `position_index` vector in the given vetcer `A`.

# Examples
```julia-repl
julia> A =  [[0.09,0.09,0.1],[0.3,0.09,0.1],[0.09]]
3-element Vector{Vector{Float64}}:
 [0.09, 0.09, 0.1]
 [0.3, 0.09, 0.1]
 [0.09]

julia> position_index =  [1,2,3]
3-element Vector{Int64}:
 1
 2
 3

julia> find_num_in_vec(A::Vector, position_index::Vector{Int64}, 0.09)
3-element Vector{Int64}:
2
1
1
```
"""
@inline function find_num_in_vec(A::Vector{Vector{T}}, position_index::Vector{Int64}, x::T) where {T}
    number_in_vec = Vector{Int64}(undef, length(position_index))
    @inbounds for i in eachindex(position_index)
        number_in_vec[i] = count(==(x),A[position_index[i]])
    end
    return number_in_vec
end

"""
    find_next_delay_vec(A::Vector{Vector{T}}, x::T)

Returns two vectors. The first is the position index vector `position_index` of vector `A`, and the second is the vector `num_in_vec` composed of the number of values in the position index corresponding to `position_index` equal to `x`.

find the minimal dt_delay in various delay channel.

# Examples
```julia-repl
julia> A =  [[0.09,0.09,0.1],[0.3,0.09,0.1],[0.09]]
3-element Vector{Vector{Float64}}:
 [0.09, 0.09, 0.1]
 [0.3, 0.09, 0.1]
 [0.09]

 julia> x=0.09
 0.09

julia> find_next_delay_vec(A, x)
([1, 2, 3], [2, 1, 1])
```
"""
function find_next_delay_vec(A::Vector{Vector{T}}, x::T) where {T}
    position_index = findall(A->isequal_in_vec(A,x),A)
    num_in_vec = find_num_in_vec(A, position_index, x)
    return position_index, num_in_vec
end


de_chan = [[0.09,0.09,0.1],[0.3,0.09,0.1],[0.09]]

find_next_delay_vec(de_chan, 0.09)

using BenchmarkTools
@benchmark find_next_delay_vec(de_chan, 0.09)
1
# @inbounds 不检查索引时的范围，加入只提升了一点点效果。
# @fastmath 使用不安全的浮点数计算以加快速度，但是计算很少，加上也没有效果。不添加。
# @inline  用于函数内联, 因为是编译器自动内联的，基本上自动内联都已经完成了，这个宏只是推荐，所以加上效果也不明显。
