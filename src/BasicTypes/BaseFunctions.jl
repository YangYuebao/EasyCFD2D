# 离散格式的总类型，
# 对流格式是离散格式的一个子抽象类型，
# 对流格式在FlowSolver的FlowTerm.jl中定义
abstract type AbstractScheme end

struct Point5 <:AbstractScheme end          # 5点格式
struct Point9 <:AbstractScheme end          # 交叉9点格式
struct Point9Orthor <:AbstractScheme end    # 正交9点区域

# getPlace用来获取不同离散格式的离散区域
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

# 根据离散格式，生成稀疏矩阵的索引
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

# 坐标转换
function k2ij(n::Int64, m::Int64, k::Int64)
    j::Int64 = 1 + floor(Int64, k / n)
    i::Int64 = k % n
    return [i, j]
end
function ij2k(n::Int64, m::Int64, i::Int64, j::Int64)
    i + (j - 1) * n
end

# 根据边界的i，j坐标求出该点处的边界条件中A，B，C分别是多少
# 边界条件：Aϕ+B(∂ϕ/∂n)=C
function toABC(n::Int64,m::Int64,i::Int64, j::Int64,alpha::Matrix{Float64},gamma::Matrix{Float64}, phi_symbol::Symbol, bounds::Vector{bound})
    if i == 1 && j == 1
        temp1 = -getproperty(bounds[1], phi_symbol)(bounds[1].span[1])
        temp4 = getproperty(bounds[4], phi_symbol)(bounds[4].span[2])
        if temp1[2] == 0
            return temp1
        end
        if temp4[2] == 0
            return temp4
        end
        return [
            temp1[1] * temp4[2] * sqrt(alpha[1, 1]) - temp4[1] * temp1[2] * sqrt(gamma[1, 1]),
            temp1[2] * temp4[2],
            temp1[3] * temp4[2] * sqrt(alpha[1, 1]) - temp4[3] * temp1[2] * sqrt(gamma[1, 1]),
        ]
    end
    if i == n && j == 1
        temp1 = -getproperty(bounds[1], phi_symbol)(bounds[1].span[2])
        temp2 = -getproperty(bounds[2], phi_symbol)(bounds[2].span[1])
        if temp1[2] == 0
            return temp1
        end
        if temp2[2] == 0
            return temp2
        end
        return [
            temp1[1] * temp2[2] * sqrt(alpha[n, 1]) + temp2[1] * temp1[2] * sqrt(gamma[n, 1]),
            temp1[2] * temp2[2],
            temp1[3] * temp2[2] * sqrt(alpha[n, 1]) + temp2[3] * temp1[2] * sqrt(gamma[n, 1]),
        ]
    end
    if i == n && j == m
        temp3 = getproperty(bounds[3], phi_symbol)(bounds[3].span[1])
        temp2 = -getproperty(bounds[2], phi_symbol)(bounds[2].span[2])
        if temp3[2] == 0
            return temp3
        end
        if temp2[2] == 0
            return temp2
        end
        return [
            temp2[1] * temp3[2] * sqrt(alpha[n, m]) - temp3[1] * temp2[2] * sqrt(gamma[n, m]),
            temp3[2] * temp2[2],
            temp2[3] * temp3[2] * sqrt(alpha[n, m]) - temp3[3] * temp2[2] * sqrt(gamma[n, m]),
        ]
    end
    if i == 1 && j == m
        temp3 = getproperty(bounds[3], phi_symbol)(bounds[3].span[2])
        temp4 = getproperty(bounds[4], phi_symbol)(bounds[4].span[1])
        if temp3[2] == 0
            return temp3
        end
        if temp4[2] == 0
            return temp4
        end
        return -[
            temp4[1] * temp3[2] * sqrt(alpha[1, m]) + temp3[1] * temp4[2] * sqrt(gamma[1, m]),
            temp3[2] * temp4[2],
            temp4[3] * temp3[2] * sqrt(alpha[1, m]) + temp3[3] * temp4[2] * sqrt(gamma[1, m]),
        ]
    end
    if i == 1
        t = (bounds[4].span[2] - bounds[4].span[1]) * (j - 1) / (m - 1)
        return getproperty(bounds[4], phi_symbol)(t)
    end
    if i == n
        t = (bounds[2].span[2] - bounds[2].span[1]) * (j - 1) / (m - 1)
        return -getproperty(bounds[2], phi_symbol)(t)
    end
    if j == 1
        t = (bounds[1].span[2] - bounds[1].span[1]) * (i - 1) / (n - 1)
        return -getproperty(bounds[1], phi_symbol)(t)
    end
    if j == m
        t = (bounds[3].span[2] - bounds[3].span[1]) * (i - 1) / (n - 1)
        return getproperty(bounds[3], phi_symbol)(t)
    end
end

# 判断一点周围的索引是否都放入系数矩阵（边界点周围有点无法放进去）
function to_val_index(n::Int64,m::Int64,i::Int64, j::Int64, scheme::AbstractScheme)
    place = getPlace(scheme)
    t = length(place)
    res = ones(Bool, t)
    for k = 1:t
        place[k] += [i, j]
        if place[k][1] < 1 || place[k][1] > n || place[k][2] < 1 || place[k][2] > m
            res[k] = false
        end
    end
    return res
end

# 根据速度场u v生成 U V
function getUV(u::Matrix{Float64}, v::Matrix{Float64}, x_u::Matrix{Float64}, x_v::Matrix{Float64}, y_u::Matrix{Float64}, y_v::Matrix{Float64})
    U = u .* y_v - v .* x_v
    V = v .* x_u - u .* y_u
    return U, V
end
