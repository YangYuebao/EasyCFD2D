abstract type AbstractScheme end

struct Point5 <:AbstractScheme end
struct Point9 <:AbstractScheme end
struct Point9Orthor <:AbstractScheme end #正交9点区域

function getPlace(::Point5)
    [[0,-1],[-1,0],[0,0],[1,0],[0,1]]
end

function getPlace(::Point9)
    place=Vector{Vector{Int64}}(undef,0)
    for j=-1:1
        for i=-1:1
            push!(place,[i,j])
        end
    end
    return place
end

function getPlace(::Point9Orthor)
    [[0,-2],[0,-1],[-2,0],[-1,0],[0,0],[1,0],[2,0],[0,1],[0,2]]
end


function index_generation(n::Int64, m::Int64, place::Vector{Vector{Int64}})
    row = Vector{Int64}(undef, 0) # 不变
    col = Vector{Int64}(undef, 0) # 不变
    for j = 1:m
        for i = 1:n
            for move in place
                temp = [i, j] + move
                if temp[1] >= 1 && temp[1] <= n && temp[2] >= 1 && temp[2] <= m
                    push!(row, ij2k(n, m, i, j))
                    push!(col, ij2k(n, m, temp[1], temp[2]))
                end
            end
        end
    end
    return row, col
end

function k2ij(n::Int64, m::Int64, k::Int64)
    j::Int64 = 1 + floor(Int64, k / n)
    i::Int64 = k % n
    return [i, j]
end

function ij2k(n::Int64, m::Int64, i::Int64, j::Int64)
    i + (j - 1) * n
end

